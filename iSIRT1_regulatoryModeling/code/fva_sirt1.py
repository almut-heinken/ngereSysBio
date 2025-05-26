from mewpy.io import Reader, Engines, read_model
from mewpy.germ.analysis import get_real_initial_state, ifva_parallel
import os
import pandas as pd
import multiprocessing
import logging


##config

def main():
    _fraction = 0.955
    fraction_text = str(_fraction).split(".")[-1]
    root_save_folder = f"..\\results\\SIRT1\\fva_experiments\\fva_{fraction_text}_NEW"

    if not os.path.exists(root_save_folder):
        os.makedirs(root_save_folder)
    # Set up logging to file
    logging.basicConfig(filename=f'{root_save_folder}\\log.txt', level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    cpu_count = os.cpu_count() - 1
    num_processes = max(1, cpu_count)
    print(f"Using {num_processes} cpu")
    print(f"ifva_parallel with fraction={_fraction}")
    logging.info("FVA ifva_parallel")
    logging.info(f"Using {num_processes} cpu")
    logging.info(f"ifva_parallel with fraction={_fraction}")
    logging.info(f"Saving results into {root_save_folder}")

    if not os.path.exists(root_save_folder):
        os.makedirs(root_save_folder)
    sbml_fname = "..\\data\\models\\sirt1_recon3d\\Recon3D_SIRT1_generic_gpr_fix_HAM_medium_sinks_closed_irreversibles.xml"
    trn_fname = "..\\data\\pypath\\grouped_sirt1_trn_manually_curated.csv"

    fluxes = pd.read_excel("..\\data\\nutrition\\NewDiets.xlsx").set_index("Unnamed: 0").to_dict()["Basic medium"]
    irreversible = pd.read_excel("..\\data\\nutrition\\irreversible.xlsx", header = None)[0].to_list()
    # loading E. coli iMC1010 model_sirt1
    recon_reader = Reader(Engines.MetabolicSBML, sbml_fname)
    trn_reader = Reader(Engines.BooleanRegulatoryCSV,
                        trn_fname, sep=',', id_col=0, rule_col=1, header=0)

    for i in [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]:
        #i /= 10.0  # Convert i to a float between 0.0 and 1.0 with 0.1 steps
        print("i = ",i)
        logging.info(f"SIRT1 expression: {i}")
        print("Reading model_sirt1...", end="")
        model_sirt1 = read_model(recon_reader, trn_reader)
        print("OK!")
        for rxn in model_sirt1.reactions.keys():
            if "sink" in rxn:
                model_sirt1.get(rxn).bounds = (0,0)
            elif rxn in irreversible:
                model_sirt1.get(rxn).bounds = (0,1000)
        for ex in model_sirt1.exchanges.keys():
            if not ex in fluxes:
                model_sirt1.get(ex).bounds = (0,1000)
            else:
                model_sirt1.get(ex).bounds = (-1*fluxes[ex],1000)
        
        model_sirt1.get("biomass_reaction").bounds = (0,1000)
        model_sirt1.objective = {"biomass_reaction":1.0}

        logging.info(f"Optimizing for objective: {model_sirt1.objective}")
        # Create real_initial_state for the current value of i
        real_initial_state = get_real_initial_state(model_sirt1, initial_state={"SIRT1": i}, strategy="mean")

        # Perform ifva
        print("Performing FVA...", end="")
        results = ifva_parallel(model_sirt1, initial_state=real_initial_state, fraction=_fraction, num_processes=num_processes, model_path=sbml_fname,trn_path=trn_fname)
        print("OK!")
        # Save the result to a CSV file
        results.to_csv(f'{root_save_folder}\\SIRT1_{i:.1f}.csv')

if __name__ == "__main__":
    multiprocessing.freeze_support()  # For Windows support
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    main()