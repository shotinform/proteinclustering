import pandas as pd

file_path = '/home/korisnik/Desktop/ip2/ip2/supplimenatry_data-1_elac054.xlsx'
sheet_names = ['Protein associated which CVD', 'Protein associated which T2D', 'Protein associated which PTA']

df1 = pd.read_excel(file_path, sheet_name=sheet_names[0])
df2 = pd.read_excel(file_path, sheet_name=sheet_names[1])
df3 = pd.read_excel(file_path, sheet_name=sheet_names[2])

genes1 = set(df1['Gene'])
genes2 = set(df2['Gene'])
genes3 = set(df3['Gene'])

common_genes = genes1.intersection(genes2).intersection(genes3)
common_genes_df = pd.DataFrame(list(common_genes), columns=['Common Genes'])
common_genes_df.to_excel('common_genes.xlsx', index=False)
