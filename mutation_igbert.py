from transformers import BertModel, BertTokenizer, BertForMaskedLM
import torch
import argparse

def predict_masked_token(H, L, mask_position):
    # Load tokenizer and model
    tokenizer = BertTokenizer.from_pretrained("Exscientia/IgBert", do_lower_case=False)
    model = BertForMaskedLM.from_pretrained("Exscientia/IgBert")

    # Combine heavy and light chains
    sequence = ' '.join(H[0]) + " [SEP] " + ' '.join(L[0])
    
    # Tokenize the sequence
    tokens = tokenizer(
        sequence, 
        add_special_tokens=True, 
        return_tensors="pt"
    )
    
    input_ids = tokens["input_ids"]
    
    # Ensure the selected position is within bounds
    if mask_position >= input_ids.shape[1]:
        raise ValueError("Selected mask position is out of bounds!")
    
    # Get the original token at the masked position (for reference)
    original_token = tokenizer.decode([input_ids[0, mask_position].item()])
    
    # Replace the token at the chosen position with [MASK]
    input_ids[0, mask_position] = tokenizer.mask_token_id
    
    # Predict the masked token
    with torch.no_grad():
        outputs = model(input_ids)
        predictions = outputs.logits
    
    # Get the most probable token ID at the masked position
    predicted_token_id = predictions[0, mask_position].argmax().item()
    predicted_token = tokenizer.decode([predicted_token_id])
    
    # Replace the masked token with the predicted token in the original sequence
    input_ids[0, mask_position] = predicted_token_id
    
    # Decode the entire sequence with the predicted token
    predicted_sequence = tokenizer.decode(input_ids[0], skip_special_tokens=False)
    
    # Remove special tokens if needed (optional)
    predicted_sequence = predicted_sequence.replace('[CLS] ', '').replace(' [SEP]', '')
    
    return predicted_sequence


parser = argparse.ArgumentParser(description="Predict masked tokens in protein sequences using IgBert.")
parser.add_argument("--H", type=str, help="Heavy chain sequence")
parser.add_argument("--L", type=str, help="Light chain sequence")
parser.add_argument("--mask_position", type=int, help="Position to mask (0-based index)")
parser.add_argument("--jobid", type=str, help="Job ID for output file")
args = parser.parse_args()

# Convert sequences to list format (assuming they're space-separated)
H_sequence = [args.H.split()]
L_sequence = [args.L.split()]

result = predict_masked_token(H_sequence, L_sequence, args.mask_position)

with open("completed_sequence_"+args.jobid+".txt", "w") as f:
    f.write(result)
    f.close()

