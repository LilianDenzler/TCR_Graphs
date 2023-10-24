import torch
import esm
import gc
import pickle
import sys
def embed(input_data,output_file, extra_dict=None):
    with open(input_data, 'rb') as f:
        input_data = pickle.load(f)
    if extra_dict:
        with open(extra_dict, 'rb') as f:
            extra_dict = pickle.load(f)
    # Load ESM-2 model
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()  # disables dropout for deterministic results

    import torch
    model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm2_t33_650M_UR50D")
    data= input_data[:10]
    for a in range(10,len(input_data)+10,10):

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
            if extra_dict:
                 txt=id+"; "+str(extra_dict[id][0])+"; "+str(key)+"\n"
            else:   
                txt=id+"; "+str(key)+"\n"
            with open(output_file, 'a') as f:
                f.write(txt)
        gc.collect()
        data= input_data[a:a+10]
        
if __name__ == "__main__":
    embed(sys.argv[1],sys.argv[2], sys.argv[3])