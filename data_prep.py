import pandas as pd
import numpy as np
import os

#LOAD SEQUENCE DATA INTO DATA FRAMES
def load_Tom_data_ab(path=None):
    if path:
        data_file = path
    else:
        base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        data_file = os.path.join(base_path, 'data','TmData_tom_ab.txt')
    abTm_df=pd.read_csv(data_file ,sep="|")
    abTm_df.columns=[i.replace(" ","") for i in abTm_df.columns]
    for i in abTm_df.columns:
        abTm_df[i].str.strip()
    tmlist=abTm_df["tm"].to_list()
    tm2list=[]
    nancounter=0
    for i in tmlist:
        if i=="nan":
            tm2list.append(np.nan)
            nancounter+=1
        else:
            try:
                tm2list.append(float(i))
            except:
                #print("is nan?", i)
                tm2list.append(np.nan)
                nancounter+=1

    abTm_df["tm"]=tm2list
    abTm_df["lightheavy"] = abTm_df["light"] + abTm_df["heavy"]
    #print(abTm_df.shape[0])
    print("number of deleted rows:",abTm_df.shape[0] - abTm_df.dropna().shape[0])
    abTm_df= abTm_df.dropna()
    print("total usable rows", abTm_df.shape[0])
    abTm_df['light'] = abTm_df['light'].apply(lambda x: x.strip())
    abTm_df['heavy'] = abTm_df['heavy'].apply(lambda x: x.strip())
    return abTm_df
    
def load_immunocore_df(path=None):
    if path:
        data_file = path
    else:
        base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        data_file = os.path.join(base_path, 'data','TCR_SEQ_TM_data.xlsx')
    Tm_df = pd.read_excel(data_file)
    Tm_df["ALPHABETA_seq"] = Tm_df["ALPHA_seq"] + Tm_df["BETA_seq"]
    #print(Tm_df.shape[0])
    print("number of deleted rows:",Tm_df.shape[0] - Tm_df.dropna().shape[0])
    Tm_df= Tm_df.dropna()
    print("total usable rows",Tm_df.shape[0])
    Tm_df=Tm_df.rename(columns={"Tm (°C)":"tm"})
    Tm_df['ALPHA_seq'] = Tm_df['ALPHA_seq'].apply(lambda x: x.strip())
    Tm_df['BETA_seq'] = Tm_df['BETA_seq'].apply(lambda x: x.strip())
    return Tm_df

def load_TCR_AB(path_ab=None, path_tcr=None):
    abTm_df=load_Tom_data_ab(path=path_ab)
    Tm_df=load_immunocore_df(path=path_tcr)
    Tm_df['Type'] = 'TCR'
    abTm_df['Type'] = 'Ab'
    concat_df=pd.concat([abTm_df,Tm_df]).reset_index(drop=True)
    return Tm_df, abTm_df, concat_df

