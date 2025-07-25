import pandas as pd
def load_gene_names(id_mapping_file):
    df = pd.read_csv(id_mapping_file, sep='\t', header=0)
    #keep the 'reviewed' uniprot to Gene Names mapping only.
    df = df[df['Reviewed']=='reviewed']
    ## keep only the first gene for each UniProt ID
    uniprot_to_gene = {p: genes.split(' ')[0] for p, genes in zip(df['Entry'], df['Gene Names'].astype(str))}
    if 'Protein names' in df.columns:
        uniprot_to_prot_names = dict(zip(df['Entry'], df['Protein names'].astype(str)))
        #node_desc = {n: {'Protein names': uniprot_to_prot_names[n]} for n in uniprot_to_prot_names}

    uniprot_to_gene['P01019-PRO_0000420660'] = uniprot_to_gene['P01019']
    uniprot_to_prot_names['P01019-PRO_0000420660'] = uniprot_to_prot_names['P01019']
    return uniprot_to_gene, uniprot_to_prot_names