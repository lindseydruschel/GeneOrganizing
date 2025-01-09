#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os

# Set working directory
os.chdir("C:/Users/druschel/Documents/Python/Gene Homolog")
print("Working directory set to:", os.getcwd())

import pandas as pd
import mygene

# Initialize MyGene client
mg = mygene.MyGeneInfo()


# In[2]:


# Define species and taxonomic ID mapping
species_map = {"Mouse": 10090, "Human": 9606, "Rat": 10116}

# Set starting and target species
starting_species = "Rat"  # Starting species
target_species = "Mouse"  # Target species

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

# Ensure all entries in 'gene_id' are strings
data["gene_id"] = data["gene_id"].astype(str)

# Preview the first few rows of the DataFrame
data.head()


# In[4]:


# Query MyGene.info for homologs and aliases
results = mg.querymany(
    data["gene_id"].tolist(),
    scopes="symbol",
    fields="homologene,alias",
    species=starting_species_id,
    verbose=True
)

# Extract initial homolog IDs and aliases for unmatched genes
def extract_homologs_and_aliases(results, target_species_id):
    """
    Extract homolog IDs for matched genes and collect aliases for unmatched genes.

    Args:
        results (list): List of query results from MyGene.
        target_species_id (int): Taxonomic ID of the target species.

    Returns:
        dict: Dictionary mapping input gene names to homolog IDs.
        dict: Dictionary mapping unmatched gene names to their aliases.
    """
    homolog_dict = {}
    alias_dict = {}
    
    for res in results:
        query = res.get("query")  # The queried gene name
        homolog = ""  # Default to blank if no match is found
        
        # Check for homologene field
        if "homologene" in res and "genes" in res["homologene"]:
            homolog_list = res["homologene"]["genes"]
            homolog = next((str(g[1]) for g in homolog_list if g[0] == target_species_id), "")
        
        # If no homolog is found, collect aliases
        if not homolog:
            aliases = res.get("alias", [])
            alias_dict[query] = aliases if isinstance(aliases, list) else [aliases]
        
        homolog_dict[query] = homolog

    return homolog_dict, alias_dict

# Process query results to extract homologs and aliases
homolog_id_dict, alias_dict = extract_homologs_and_aliases(results, target_species_id)

# Debug print: Check the size and sample of alias_dict
print(f"Alias dictionary created with {len(alias_dict)} entries.")
print(f"Sample aliases: {list(alias_dict.items())[:5]}")



# In[5]:


import concurrent.futures

# Initialize homolog_dict_aliases (this will contain alias homologs)
homolog_dict_aliases = {}

# Step 2: Query aliases for homologs in parallel, with print statements when a homolog is found
def find_homologs_via_aliases_parallel(alias_dict, target_species_id, max_aliases=3):
    """
    Query MyGene.info for homologs using aliases of unmatched genes in parallel,
    with a limit on the number of aliases queried and print when homolog is found.

    Args:
        alias_dict (dict): Dictionary mapping unmatched gene names to their aliases.
        target_species_id (int): Taxonomic ID of the target species.
        max_aliases (int): Maximum number of aliases to query per gene.

    Returns:
        dict: Dictionary mapping unmatched gene names to homologs found via aliases.
    """
    alias_homolog_dict = {}

    # Function to query for a single alias
    def query_alias_for_gene(gene, aliases):
        for alias in aliases[:max_aliases]:  # Limit the number of aliases queried
            alias_results = mg.querymany(
                [alias],
                scopes="symbol",
                fields="homologene",
                species=target_species_id,
                verbose=False
            )
            
            for res in alias_results:
                if "homologene" in res and "genes" in res["homologene"]:
                    homolog_list = res["homologene"]["genes"]
                    for entry in homolog_list:
                        species_id, homolog_id = entry
                        if species_id == target_species_id:
                            print(f"Gene: {gene} found homolog {homolog_id} for alias: {alias}")
                            return gene, str(homolog_id)
        return gene, ""  # Return blank if no homolog found

    # Parallelize the queries for all genes
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(query_alias_for_gene, gene, aliases) for gene, aliases in alias_dict.items()]
        
        for future in concurrent.futures.as_completed(futures):
            gene, homolog = future.result()
            alias_homolog_dict[gene] = homolog

    return alias_homolog_dict

# Query aliases for unmatched genes to find homologs, with prints for homologs found
alias_homolog_dict = find_homologs_via_aliases_parallel(alias_dict, target_species_id)

# Debug print: Check the size and sample of alias_homolog_dict
print(f"Alias homolog dictionary created with {len(alias_homolog_dict)} entries.")
print(f"Sample entries: {list(alias_homolog_dict.items())[:5]}")


# In[16]:


# Step 1: Fetch Gene Names for the homolog IDs from the original and alias homologs
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


# Fetch gene names for original homologs (homolog_id_dict)
homolog_name_dict = fetch_gene_names(homolog_id_dict)

# Fetch gene names for alias homologs (alias_homolog_dict)
alias_name_dict = fetch_gene_names(alias_homolog_dict)

# Step 2: Create a new dictionary mapping original `gene_id` to the final homolog name
final_homolog_dict = {}

# Add original homolog names first
for gene, homolog_id in homolog_id_dict.items():
    if homolog_id in homolog_name_dict:  # Check if a gene name exists for this homolog ID
        final_homolog_dict[gene] = homolog_name_dict[homolog_id]

# Add alias homolog names only if no original homolog exists
for gene, alias_homolog_id in alias_homolog_dict.items():
    if gene not in final_homolog_dict and alias_homolog_id in alias_name_dict:
        final_homolog_dict[gene] = alias_name_dict[alias_homolog_id]

# Debug: Check the final homolog dictionary
print("\nFinal homolog dictionary mapping gene IDs to names:")
print(f"Sample from final_homolog_dict: {list(final_homolog_dict.items())[:5]}")

# Step 3: Update DataFrame with the final homolog names (with case matching)

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
    match_case(gene, final_homolog_dict.get(gene, ""))
    for gene in data["gene_id"]
]

# Add aliases to the DataFrame (with case matching)
data["Aliases"] = [
    ", ".join([match_case(gene, alias) for alias in alias_dict.get(gene, []) if isinstance(alias, str)])
    for gene in data["gene_id"]
]

# Preview the updated DataFrame
print("\nUpdated DataFrame preview:")
print(data.head())

# Calculate match rate
total_genes = len(data)
matched_genes = len([gene for gene in data["gene_id"] if final_homolog_dict.get(gene, "") != ""])
match_rate = (matched_genes / total_genes) * 100
print(f"\nMatch rate: {matched_genes}/{total_genes} genes matched ({match_rate:.2f}%)")






# In[18]:


# Save all results to a CSV file
output_csv_all = "output_genes_with_aliases_all.csv"
data.to_csv(output_csv_all, index=False)
print(f"All results saved to: {output_csv_all}")

# Save matched results to a separate CSV file
matched_data = data[data["New_Homolog"] != ""]
output_csv_matched = "output_genes_with_aliases_matched.csv"
matched_data.to_csv(output_csv_matched, index=False)
print(f"Matched results saved to: {output_csv_matched}")

# Calculate match rate
total_genes = len(data)
matched_genes = len(matched_data)
match_rate = (matched_genes / total_genes) * 100
print(f"Match rate: {matched_genes}/{total_genes} genes matched ({match_rate:.2f}%)")


# In[ ]:





# In[ ]:




