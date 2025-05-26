import cobra
import multiprocessing
import pandas as pd
import re
import os
import csv

from cobra.io import read_sbml_model, write_sbml_model
from corda import reaction_confidence
from pyfastcore import Fastcore
from copy import deepcopy


config = cobra.Configuration()
config.solver = 'cplex'



# Define the filename for CSV import
csv_filename = "..\\results\\GTEx\\data\\core_reactions.csv"

# Initialize an empty list to store the loaded data
loaded_core_reactions = []

# Load the data from CSV
with open(csv_filename, mode='r') as file:
    reader = csv.reader(file)
    # Skip the header row if present
    next(reader, None)
    # Read each row and append the data to the list
    for row in reader:
        loaded_core_reactions.append(row[0])

print("Loaded core reactions from CSV:", loaded_core_reactions)

dataset = '..\\data\\GTEx\\RNASeq_log2_normalized.parquet'
sbml_fname = "..\\data\\models\\sirt1_recon3d\\Recon3D_SIRT1_generic_gpr_fix_HAM_medium_sinks_closed_irreversibles.xml"

#Read FVA results for raw Recon3D with fraction of 0.955
fva_file = "..\\results\\SIRT1\\fva_experiments\\fva_955_NEW\\SIRT1_0.0.csv"
fva_00 = pd.read_csv(fva_file).set_index("Unnamed: 0")

output_folder = "..\\results\\GTEx\\fastcore_models_HAM_medium_sinks_closed_irreversibles_core\\"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

quantile_th = [0.1,0.5,0.75]

gene_expression_df = pd.read_parquet(dataset)

print("Reading model...", end="")
model = read_sbml_model(sbml_fname)

# for rxn in fva_00.index.to_list():
#     if model.reactions.get_by_id(rxn):
#         _minimum = fva_00.loc[rxn].minimum
#         _maximum = fva_00.loc[rxn].maximum
#         model.reactions.get_by_id(rxn).bounds = (_minimum, _maximum)

# model.reactions.get_by_id("biomass_reaction").bounds = (0,1000)
print("OK!")

global model_genes
model_genes = [g.id for g in model.genes]

def create_rxn_confidence(gene_confidence_dict):
    rxn_conf_dict = {}
    for r in model.reactions:
        if re.match('biomass_reaction', r.id) or r.id in loaded_core_reactions:
            rxn_conf_dict[r.id] = 3
        else:
            if len(r.genes) > 0:
                rxn_conf_dict[r.id] =  reaction_confidence(r, gene_confidence_dict)
            else:
                rxn_conf_dict[r.id] =  0
    return rxn_conf_dict


def gene_confidence(gene_dict,col):
    gene_confidence_dict = {}
    confidence_zero_genes = set(model_genes) - set(gene_dict.keys())
    
    for gene in confidence_zero_genes:
        gene_confidence_dict[gene] = 0
        
    for k in model_genes:
        try:
            v = gene_dict[k]
        except KeyError:
            continue
        
        gene_expr = gene_expression_df[col]
        q_1 = gene_expr.quantile(quantile_th[0])
        q_2 = gene_expr.quantile(quantile_th[1])
        q_4 = gene_expr.quantile(quantile_th[2])
        if not v == 'nan':
            if v <= q_1:
                gene_confidence_dict[k] = -1
            elif v <= q_2:
                gene_confidence_dict[k] = 1
            elif v < q_4:
                gene_confidence_dict[k] = 2
            elif v >= q_4:
                gene_confidence_dict[k] = 3
        else:
            #If 'nan' then assign confidence == 0
            gene_confidence_dict[k] = 0
    
    return gene_confidence_dict



def run_fastcore(gene_dict_wt):
    column_name = gene_dict_wt['column_name']
    print(f"Computing column: {column_name}", flush=True)
    gene_confidence_wt_dict = gene_confidence(gene_dict_wt['values'],column_name)
    reactions_confidence_wt_dict = create_rxn_confidence(gene_confidence_wt_dict)
    core_reactions = [r for r,v in reactions_confidence_wt_dict.items() if v > 2]
    print("core reactions:",len(core_reactions))
    # Setting the penalty of exchange fluxes to 0
    penalties = {}
    for r in model.exchanges:
        penalties[r.id] = 0

    # Creating a fastcore solver instnace
    model_wt = deepcopy(model)
    fc_builder = Fastcore(model_wt, core_reactions,
                          penalties=penalties,
                          default_penalty=10,
                          debug_mode=True)

    # Rnunning fastcore
    fc_builder.fast_core()
    
    consistent_subnetwork = fc_builder.consistent_subnetwork
    print(f"Consistent sub-network lenght: {len(consistent_subnetwork)}")
    
    cs_model = fc_builder.build_context_specific_model()

    print(f"Saving model as Recon3D_{column_name}.xml...")
    write_sbml_model(cs_model,output_folder+f"Recon3D_{column_name}.xml")



if __name__ == "__main__":
    # Create a list of dictionaries, each containing column name and values
    gene_dicts = [{'column_name': col, 'values': gene_expression_df.to_dict()[col]} for col in gene_expression_df.columns]

    # Create a multiprocessing pool
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()-2)  # Use all available CPU cores

    # Map the run_fastcore function to process each column in parallel
    pool.map(run_fastcore, gene_dicts)

    # Close the pool to release resources
    pool.close()
    pool.join()