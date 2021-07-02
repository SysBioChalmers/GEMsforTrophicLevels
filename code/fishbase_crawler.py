import pandas as pd

table = {}
table_tmp = pd.read_html('https://www.fishbase.de/trophic/T.html')
table_tmp = [x for x in table_tmp if x.shape[1] == 6]
table_tmp = pd.concat(table_tmp)
table_tmp = table_tmp.iloc[:,1:].reset_index(drop = True)
table_tmp = table_tmp.rename(columns = table_tmp.iloc[0]).iloc[1:]
table[0] = table_tmp

for i in range(2, 74):
    table_tmp = pd.read_html(f"https://www.fishbase.de/trophic/TPage{i}.html")
    table_tmp = [x for x in table_tmp if x.shape[1] == 6]
    table_tmp = pd.concat(table_tmp)
    table_tmp = table_tmp.iloc[:,1:].reset_index(drop = True)
    table_tmp = table_tmp.rename(columns = table_tmp.iloc[0]).iloc[1:]
    table[i] = table_tmp

table_full = pd.concat(table)
table_full.reset_index(drop = True, inplace = True)
table_full.to_csv('Fish trophic levels')