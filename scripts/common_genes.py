import pandas as pd

# skripta koja pronalazi zajedničke gene za sve tri bolesti
directory_path = '../files/disease_genes/'
file_names = ['genes_PTA.xlsx', 'genes_T2D.xlsx', 'genes_CVD.xlsx']

df1 = pd.read_excel(directory_path+file_names[0])
df2 = pd.read_excel(directory_path+file_names[1])
df3 = pd.read_excel(directory_path+file_names[2])

genes1 = set(df1['Gene'])
genes2 = set(df2['Gene'])
genes3 = set(df3['Gene'])

output_file = '../files/common_genes.xlsx'
common_genes = genes1.intersection(genes2).intersection(genes3)
common_genes_df = pd.DataFrame(list(common_genes), columns=['Common Genes'])
common_genes_df.to_excel(output_file, index=False)
