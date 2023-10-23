# read csv
import pandas as pd
wt_group_df=pd.read_csv("/home/lilian/TCR_Graphs/Novozyme_data/AF_WT/train_wildtype_groups.csv")
#count how many rows have same value for ph
df=wt_group_df["pH"].value_counts()
#keep only rows with ph=7
#wt_group_df=wt_group_df[wt_group_df["pH"]==7]
wt_group_df=wt_group_df[wt_group_df["pH"]>=6.5]
wt_group_df=wt_group_df[wt_group_df["pH"]<=7.5]
wt_group_df.columns

import torch
model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm2_t33_650M_UR50D")

df=wt_group_df[["seq_id","protein_sequence"]]
df.set_index('seq_id')
df=df[["protein_sequence"]]   
dict=df.T.to_dict("list")

dict = { str(k): v for k, v in dict.items() }
dict = { k: v[0] for k, v in dict.items() }
protein_seq_list=list(dict.items())

import torch
import esm
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results

def encode(data_list, model, alphabet, batch_converter):
    # Prepare data (first 2 sequences from ESMStructuralSplitDataset superfamily / 4)
    data = protein_seq_list
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
    return sequence_representations

# Look at the unsupervised self-attention map contact predictions
"""import matplotlib.pyplot as plt
for (_, seq), tokens_len, attention_contacts in zip(data, batch_lens, results["contacts"]):
    plt.matshow(attention_contacts[: tokens_len, : tokens_len])
    plt.title(seq)
    plt.show()"""
    
#write list to file
import csv
for i in range(0,len(protein_seq_list),1):
    try:
        data=protein_seq_list[i:i+1]
        print(data)
        sequence_representations=encode(data, model, alphabet, batch_converter)
        with open('/home/lilian/TCR_Graphs/Novozyme_data/AF_WT/seq_representations_ph6575.csv', 'w') as f:
            for key in sequence_representations:
                f.write("%s\n" % key)
    except:
        print("error"  )
        data=protein_seq_list[i:]
        sequence_representations=encode(data, model, alphabet, batch_converter)
        with open('/home/lilian/TCR_Graphs/Novozyme_data/AF_WT/seq_representations_ph6575.csv', 'w') as f:
            for key in sequence_representations:
                f.write("%s\n" % key)
