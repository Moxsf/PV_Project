path_in = "/h/xutong/project_RV/21_nonSigGWAS/08_INDEL/code/"
path_out = "/h/xutong/project_RV/21_nonSigGWAS/08_INDEL/code/"

filenm = "cc"
file_in = path_in + filenm
file_out = path_out + "bb"

Data_in = pd.read_csv(file_in, sep="\t", header=0)

ids = Data_in["class"].unique()
Data_out = pd.DataFrame()

ids = Data_in["class"].unique()
Data_out = pd.DataFrame()

for idnm in ids:
    temp_data = Data_in[Data_in["class"].isin([idnm])]
    temp_max = pd.DataFrame(temp_data.max()).T
    Data_out = pd.concat([Data_out, temp_max])

Data_out["chrom"] = Data_out["class"].map(lambda x:x.split("_")[0])
Data_out["pos"] = Data_out["class"].map(lambda x:x.split("_")[1])
Data_out["ref"] = Data_out["class"].map(lambda x:x.split("_")[2])
Data_out["alt"] = Data_out["class"].map(lambda x:x.split("_")[3])


for RowNum in range(Data_out.shape[0]):
    if len(Data_out.iloc[RowNum, ]["ref"]) > len(Data_out.iloc[RowNum, ]["alt"]):
        Data_out.iloc[RowNum, ]["type"] = "DEL"
        Data_out.iloc[RowNum, ]["length"] = abs(len(Data_out.iloc[RowNum, ]["ref"]) - len(Data_out.iloc[RowNum, ]["alt"]))
    elif len(Data_out.iloc[RowNum, ]["ref"]) < len(Data_out.iloc[RowNum, ]["alt"]):
        Data_out.iloc[RowNum, ]["type"] = "INS"
        Data_out.iloc[RowNum, ]["length"] = abs(len(Data_out.iloc[RowNum, ]["ref"]) - len(Data_out.iloc[RowNum, ]["alt"]))
    else:
        pass

Data_out.to_csv(file_out, sep="\t", index=False)
