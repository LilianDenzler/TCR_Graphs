#Read in TCR sequences and antibody sequences using data_prep.py
#Run ANARCI on each sequence to get Fasta files of each correctly numbered sequence

import data_prep
import numpy as np
import pandas as pd
import subprocess
import os
from natsort import natsorted

def anarci_number(full_df,sequence, seqid, tm,  path="./anarcinumbered_imgt"):
    #if folder doesn't exist, create it
    if not os.path.exists(path):
        os.makedirs(path)

    subprocess.run(["ANARCI", "-i", sequence, "--scheme", "imgt", "--restrict", "tr", "--csv", "-o",os.path.join(path,seqid+".csv")], capture_output=False)
    try:
        Adf=pd.read_csv(os.path.join(path,seqid+".csv_A.csv"))
        Bdf=pd.read_csv(os.path.join(path,seqid+".csv_B.csv"))
    except:
        return full_df
    Adf["tm"]=[tm]
    Bdf["tm"]=[tm]
    df=pd.concat([Adf, Bdf])
    if full_df.empty:
        full_df=df
    else:
        full_df=pd.concat([full_df,df])
    return full_df
    

def numberdf(Tm_df, path="./full_df.csv"):
    full_df=pd.DataFrame()
    Tm_df["ALPHABETA_seq"]=Tm_df["ALPHABETA_seq"].apply(lambda x: x.replace(" ",""))
    Tm_df["ALPHABETA_seq"]=Tm_df["ALPHABETA_seq"].apply(lambda x: x.replace(".",""))
    Tm_df["ALPHABETA_seq"]=Tm_df["ALPHABETA_seq"].apply(lambda x: x.replace("-",""))
    for (i, id,tm) in zip(Tm_df["ALPHABETA_seq"], Tm_df["ID"], Tm_df["tm"]):
        full_df=anarci_number(full_df,i, id,tm)
    full_df.to_csv(path)
    return path

def get_dfs(full_df_path="./full_df.csv"):
    df=pd.read_csv(full_df_path)
    df = df.fillna("-")
    df.rename(columns={"Unnamed: 0": "tm"}, inplace = True)
    df.drop(df.columns[df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
    """for col in df.columns:
        if col in ["tm","Id", "sequence", "tm", "species", "chain_type", "v_gene", "j_gene", "v_identity", "j_identity", "hmm_species", "identity_species", "seqstart_index", "seqend_index", "e-value", "score", "domain_no"]:
            continue
        else:
            df[col] = df[col].apply(lambda x: aa2int[x])"""

    encoded_df = df
    Adf_encoded = df[df['chain_type'] == "A"]
    Bdf_encoded = df[df['chain_type'] == "B"]
    # Split data into training and testing sets
    Adf_tm= Adf_encoded.drop(["chain_type", "v_gene", "j_gene", "v_identity", "j_identity", "hmm_species", "identity_species", "seqstart_index", "seqend_index", "e-value", "score", "domain_no"], axis=1)
    Adf= Adf_tm.drop(["tm","Id"], axis=1)
    Bdf_tm= Bdf_encoded.drop(["chain_type", "v_gene", "j_gene", "v_identity", "j_identity", "hmm_species", "identity_species", "seqstart_index", "seqend_index", "e-value", "score", "domain_no"], axis=1)
    Bdf= Bdf_tm.drop(["tm","Id"], axis=1)
    return Adf, Adf_tm, Bdf, Bdf_tm

    
def df_to_fasta(df, df_tmId, path="./fasta"):
    if not os.path.exists(path):
        os.makedirs(path)
    if not os.path.exists(path+"_nomissing"):
        os.makedirs(path+"_nomissing")
    for i in list(df.index.values):
        i=int(i)
        row=df.loc[i, :].values.flatten().tolist()
        with open(os.path.join(path,"seq_"+str(i)+".fasta"), "w") as f:
            f.write(">"+"seq_"+str(i)+"\n")
            f.write("".join(str(x) for x in row))
        with open(os.path.join(path+"_nomissing","seq_"+str(i)+".fasta"), "w") as f:
            f.write(">"+"seq_"+str(i)+"\n")
            f.write("".join(str(x) for x in row if x != "-"))
    
def run():
    #Tm_df=data_prep.load_immunocore_df(path='./TCR_SEQ_TM_data.xlsx')
    #numberdf(Tm_df, path="./full_df.csv")
    Adf, Adf_tmId, Bdf, Bdf_tmId = get_dfs(full_df_path="./full_df.csv")
    Adf = Adf.reindex(natsorted(Adf.columns), axis=1)
    Bdf = Bdf.reindex(natsorted(Bdf.columns), axis=1)
    print(Adf)
    print(Bdf)
    df_to_fasta(Adf, Adf_tmId, path="./fastaA")
    df_to_fasta(Bdf, Bdf_tmId, path="./fastaB")
    
if __name__ == "__main__":
    run()