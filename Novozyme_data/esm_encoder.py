import torch
import esm
import gc
# Load ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results


# read csv
import pandas as pd
wt_group_df=pd.read_csv("/home/lilian/TCR_Graphs/Novozyme_data/AF_WT/train_wildtype_groups_ph_6.5_7.5.csv", index_col=None, header=0)
#count how many rows have same value for ph
df=wt_group_df["pH"].value_counts()
#keep only rows with ph=7
#wt_group_df=wt_group_df[wt_group_df["pH"]==7]
wt_group_df=wt_group_df[wt_group_df["pH"]>=6.5]
wt_group_df=wt_group_df[wt_group_df["pH"]<=7.5]
wt_group_df.columns

import torch
model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm2_t33_650M_UR50D")

df=wt_group_df[["seq_id","protein_sequence", "tm"]]
df=df.set_index('seq_id')
seq_df=df[["protein_sequence"]]
tm_df=df[["tm"]]   

tm_dict=tm_df.T.to_dict("list")  
tm_dict = { str(k): v for k, v in tm_dict.items() }

dict=seq_df.T.to_dict("list")
dict = { str(k): v for k, v in dict.items() }
dict = { k: v[0] for k, v in dict.items() }
protein_seq_list=list(dict.items())
seq_id_list=list(dict.keys())

data= protein_seq_list[:10]
for a in range(10,len(protein_seq_list),10):

    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

    # Extract per-residue representations (on CPU)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=True)
    token_representations = results["representations"][33]

    # Generate per-sequence representations via averaging
    # NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
    sequence_representations = []
    for i, tokens_len in enumerate(batch_lens):
        sequence_representations.append(token_representations[i, 1 : tokens_len - 1].mean(0))

    
    for id, key in zip(batch_labels,sequence_representations):
        print(id)
        key=key.tolist()
        txt=id+"; "+str(tm_dict[id][0])+"; "+str(key)+"\n"
        with open('/home/lilian/TCR_Graphs/Novozyme_data/AF_WT/seq_representations_ph6575.csv', 'a') as f:
            f.write(txt)
    gc.collect()
    data= protein_seq_list[a:a+10]