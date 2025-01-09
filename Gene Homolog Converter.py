#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import pandas as pd
import mygene

# Initialize MyGene client
mg = mygene.MyGeneInfo()

# Set working directory
os.chdir("C:/Users/druschel/Documents/Python/Gene Homolog")
print("Working directory set to:", os.getcwd())



# In[2]:


# Define species and taxonomic ID mapping
species_map = {"Mouse": 10090, "Human": 9606, "Rat": 10116}

# Set starting and target species
starting_species = "Rat"  # Starting species (e.g., "Rat")
target_species = "Mouse"  # Target species (e.g., "Mouse")

# Get taxonomic IDs for starting and target species
starting_species_id = species_map[starting_species]
target_species_id = species_map[target_species]

print(f"Starting Species: {starting_species} (ID: {starting_species_id})")
print(f"Target Species: {target_species} (ID: {target_species_id})")


# In[3]:


# Load the input Excel file
file_path = "Yulia all 17k important FPKM.xlsx"  # Replace with your actual file name
data = pd.read_excel(file_path)

# Ensure the file contains a "gene_id" column
if "gene_id" not in data.columns:
    raise ValueError("The input file must contain a 'gene_id' column.")

# Preview the first few rows of the DataFrame
data.head()


# In[4]:


# Extract genes from the input data
genes = data["gene_id"].tolist()

# Query homolog data from MyGene.info
results = mg.querymany(genes, scopes="symbol", fields="homologene", species=starting_species_id, verbose=True)

# Preview the first few query results
results[:5]


# In[5]:


# Step 1: Extract homolog gene IDs - takes the highest confidence homolog but does NOT filter based on confidence
def process_homolog_ids_highest_score(results, target_species_id):
    """
    Process MyGene query results to extract homolog IDs for a target species,
    selecting the homolog with the highest score if multiple are found.

    Args:
        results (list): List of query results from MyGene.
        target_species_id (int): Taxonomic ID of the target species.

    Returns:
        dict: Dictionary mapping input gene names to homolog IDs.
    """
    result_dict = {}
    
    for res in results:
        query = res.get("query")  # The queried gene name
        best_homolog = ""  # Default to blank if no match is found
        best_score = float('-inf')  # Track the highest score
        
        # Check if "homologene" exists and contains valid data
        if "homologene" in res and "genes" in res["homologene"]:
            homolog_list = res["homologene"]["genes"]
            
            for entry in homolog_list:
                species_id, homolog_id = entry
                # Check if the entry matches the target species
                if species_id == target_species_id:
                    # Get the `_score` for the result (default to 0 if missing)
                    score = res.get("_score", 0)
                    # Update the best homolog if this one has a higher score
                    if score > best_score:
                        best_score = score
                        best_homolog = str(homolog_id)
        
        # Add the best match to the result dictionary
        result_dict[query] = best_homolog

    return result_dict

# Process query results to extract homolog IDs with the highest confidence
homolog_id_dict = process_homolog_ids_highest_score(results, target_species_id)



# In[6]:


# Step 2: Query MyGene.info for gene names
def fetch_gene_names(homolog_ids):
    """
    Fetch gene names (symbols) for homolog IDs using MyGene.info.

    Args:
        homolog_ids (dict): Dictionary mapping gene names to homolog IDs.
    
    Returns:
        dict: Dictionary mapping homolog IDs to gene names.
    """
    gene_name_dict = {}
    # Filter out "No Match" and empty homolog IDs
    valid_ids = [homolog_id for homolog_id in homolog_ids.values() if homolog_id != "No Match"]
    
    # Query MyGene.info with homolog IDs
    homolog_name_results = mg.querymany(valid_ids, scopes="entrezgene", fields="symbol", species=target_species_id, verbose=True)
    
    for res in homolog_name_results:
        gene_id = res.get("query")
        gene_name = res.get("symbol", "No Name Found")  # Default to "No Name Found" if symbol is missing
        gene_name_dict[gene_id] = gene_name

    return gene_name_dict

# Fetch gene names for homolog IDs
homolog_name_dict = fetch_gene_names(homolog_id_dict)


# In[7]:


# Map homolog names to input data
data["New_Homolog"] = [homolog_name_dict.get(homolog_id_dict.get(gene, ""), "No Match") for gene in data["gene_id"]]

# Preview the updated DataFrame
data.head()


# In[8]:


# Ensure all entries in the 'gene_id' column are strings
data["gene_id"] = data["gene_id"].astype(str)

# Function to match the case style of the input
def match_case(input_name, output_name):
    if input_name.isupper():  # If the input is all uppercase
        return output_name.upper()
    elif input_name.istitle():  # If the input is title case
        return output_name.title()
    else:  # Default to lowercase if input has mixed or lowercase
        return output_name.lower()

# Map homolog names to input data while preserving case
data["New_Homolog"] = [
    match_case(gene, homolog_name_dict.get(homolog_id_dict.get(gene, ""), ""))
    for gene in data["gene_id"]
]

# Preview the updated DataFrame
data.head()



# In[9]:


## Save the updated data to an Excel file
#output_excel_file = "output_genes_with_homologs.xlsx"  # Replace with your desired file name
#data.to_excel(output_excel_file, index=False)
#print(f"Processed Excel file saved to: {output_excel_file}")


# Save all results (including blanks) to a CSV file
output_csv_all = "output_genes_all.csv"
data.to_csv(output_csv_all, index=False)
print(f"All results saved to: {output_csv_all}")

# Filter to keep only matched genes (non-blank New_Homolog)
matched_data = data[data["New_Homolog"] != ""]

# Save matched results to a separate CSV file
output_csv_matched = "output_genes_matched.csv"
matched_data.to_csv(output_csv_matched, index=False)
print(f"Matched results saved to: {output_csv_matched}")

# Calculate and print match rate
total_genes = len(data)
matched_genes = len(matched_data)
match_rate = (matched_genes / total_genes) * 100
print(f"Match rate: {matched_genes}/{total_genes} genes matched ({match_rate:.2f}%)")


# In[ ]:




