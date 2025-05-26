from mewpy.io import Reader, Engines, read_model
from mewpy.germ.analysis import get_real_initial_state, ifva_parallel
import os
import pandas as pd
import multiprocessing
import logging


##config
_fraction = 0.955

def main():
    root_save_folder = f"..\\results\\caco2\\fva_simulations"
    if not os.path.exists(root_save_folder):
        os.makedirs(root_save_folder)
    # Set up logging to file
    logging.basicConfig(filename=f'{root_save_folder}\\log.txt', level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    cpu_count = os.cpu_count() - 1
    num_processes = max(1, cpu_count)
    print(f"Using {num_processes} cpu")
    
    print(f"ifva_parallel with fraction={_fraction}")
    logging.info(f"Using {num_processes} cpu")
    logging.info(f"ifva_parallel with fraction={_fraction}")

    file_path = "..\\data\\caco2\\DMEM.xlsx"
    df = pd.read_excel(file_path, header=2).set_index("Reaction ID")
    df.columns = ["description", "uptake"]
    medium = df.to_dict()["uptake"]


    for butyrate_flux in [0.0, 0.005, 0.01, 0.03, 0.05]:
        logging.info(f"butyrate_flux: {butyrate_flux}")
        print("butyrate_flux:",butyrate_flux)
        save_folder = root_save_folder
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        sbml_fname = "..\\data\\models\\sirt1_recon3d\\Recon3D_SIRT1_generic_gpr_fix.xml"
        trn_fname = "..\\data\\pypath\\grouped_sirt1_trn_manually_curated.csv"

        # loading E. coli iMC1010 model
        recon_reader = Reader(Engines.MetabolicSBML, sbml_fname)
        trn_reader = Reader(Engines.BooleanRegulatoryCSV,
                            trn_fname, sep=',', id_col=0, rule_col=1, header=0)
        
        print("Reading model_sirt1...", end="")
        model_sirt1 = read_model(recon_reader, trn_reader)
        print("OK!")

        print("Setting flux boundaries...", end="")
        for ex in model_sirt1.exchanges:
            exchange = model_sirt1.get(ex)
            if ex in medium:
                exchange.bounds = (-1*medium[ex], 1000)
            else:
                exchange.bounds = (0, 1000)
        print("OK!")

        print("Setting sink boundaries to (0,0)...", end="")
        for rxn in model_sirt1.reactions.keys():
            if "sink" in rxn:
                model_sirt1.get(rxn).bounds = (0,0)
        print("OK!")

        #i /= 10.0  # Convert i to a float between 0.0 and 1.0 with 0.1 steps
        print(f"Setting EX_but[e] boundaries to {(-1*butyrate_flux, 1000)}...")
        logging.info(f"butyrate_flux: {butyrate_flux}")
        model_sirt1.get("EX_but[e]").bounds = (-1*butyrate_flux, 1000)
        sirt1_expression = max(0, -6.116 * butyrate_flux + 0.945)
        print(f"sirt1_expression to {sirt1_expression}...")
        logging.info(f"SIRT1 expression: {sirt1_expression}")
        
        model_sirt1.objective = {"biomass_reaction":1.0}
        logging.info(f"Optimizing for objective: {model_sirt1.objective}")
        # Create real_initial_state for the current value of i
        real_initial_state = get_real_initial_state(model_sirt1, initial_state={"SIRT1": sirt1_expression}, strategy="mean")

        # Perform ifva
        print("Performing FVA...", end="")
        results = ifva_parallel(model_sirt1, initial_state=real_initial_state, fraction=_fraction, num_processes=num_processes, model_path=sbml_fname,trn_path=trn_fname, medium=medium)
        print("OK!")
        # Save the result to a CSV file
        results.to_csv(f'{save_folder}\\butyrate_{butyrate_flux}_SIRT1_{sirt1_expression:.1f}.csv')

if __name__ == "__main__":
    multiprocessing.freeze_support()  # For Windows support
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    main()