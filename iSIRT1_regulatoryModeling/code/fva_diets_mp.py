from mewpy.io import Reader, Engines, read_model
from mewpy.germ.analysis import get_real_initial_state, ifva_parallel
import os
import multiprocessing
import logging


##config

def main():
    for _fraction in [0.955,0.99]:
        fraction_text = str(_fraction).split(".")[-1]
        root_save_folder = f"..\\results\\SIRT1\\fva_experiments\\fva_diets_{fraction_text}"
        if not os.path.exists(root_save_folder):
            os.makedirs(root_save_folder)
        # Set up logging to file
        logging.basicConfig(filename=f'{root_save_folder}\\log.txt', level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

        #diets = ["High-protein","Western", "High-carbohydrate", "High-fat", "Ketogenic"]
        diets = ["Ketogenic", "High-protein","Western", "High-carbohydrate", "High-fat"]
        #diets = ["High-protein","Western", "High-carbohydrate", "High-fat"]
        cpu_count = os.cpu_count()
        num_processes = max(1, cpu_count)
        print(f"Using {num_processes} cpu")
        
        print(f"ifva_parallel with fraction={_fraction}")
        logging.info(f"Using {num_processes} cpu")
        logging.info(f"ifva_parallel with fraction={_fraction}")

        for diet in diets:
            logging.info(f"Diet: {diet}")
            print("Diet:",diet)
            save_folder = f"{root_save_folder}\\{diet}"
            if not os.path.exists(save_folder):
                os.makedirs(save_folder)
            sbml_fname = f"..\\results\\SIRT1\\diet_models\\Recon3D_{diet}_diet.xml"
            trn_fname = "..\\data\\pypath\\grouped_sirt1_trn_manually_curated.csv"

            # loading E. coli iMC1010 model
            recon_reader = Reader(Engines.MetabolicSBML, sbml_fname)
            trn_reader = Reader(Engines.BooleanRegulatoryCSV,
                                trn_fname, sep=',', id_col=0, rule_col=1, header=0)

            for i in [0.0,0.5,1.0]:
                #i /= 10.0  # Convert i to a float between 0.0 and 1.0 with 0.1 steps
                print("i = ",i)
                logging.info(f"SIRT1 expression: {i}")
                print("Reading model_sirt1...", end="")
                model_sirt1 = read_model(recon_reader, trn_reader)
                print("OK!")
                for rxn in model_sirt1.reactions.keys():
                    if "sink" in rxn:
                        model_sirt1.get(rxn).bounds = (0,0)
                model_sirt1.objective = {"biomass_reaction":1.0}
                logging.info(f"Optimizing for objective: {model_sirt1.objective}")
                # Create real_initial_state for the current value of i
                real_initial_state = get_real_initial_state(model_sirt1, initial_state={"SIRT1": i}, strategy="mean")

                # Perform ifva
                print("Performing FVA...", end="")
                results = ifva_parallel(model_sirt1, initial_state=real_initial_state, fraction=_fraction, num_processes=num_processes, model_path=sbml_fname,trn_path=trn_fname)
                print("OK!")
                # Save the result to a CSV file
                results.to_csv(f'{save_folder}\\{diet}_SIRT1_{i:.1f}.csv')

if __name__ == "__main__":
    multiprocessing.freeze_support()  # For Windows support
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    main()