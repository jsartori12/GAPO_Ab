#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 17:05:04 2025

@author: joao
"""

from transformers import BertModel, BertTokenizer, BertForMaskedLM
import torch

# Load tokenizer and model
tokenizer = BertTokenizer.from_pretrained("Exscientia/IgBert", do_lower_case=False)
model = BertForMaskedLM.from_pretrained("Exscientia/IgBert")

# Example heavy and light chain sequences
sequences_heavy = ["VQLAQSGSELRKPGASVKVSCDTSGHSFTSNAIHWVRQAPGQGLEWMGWINTDTGTPTYAQGFTGRFVFSLDTSARTAYLQISSLKADDTAVFYCARERDYSDYFFDYWGQGTLVTVSS"]
sequences_light = ["EVVMTQSPASLSVSPGERATLSCRARASLGISTDLAWYQQRPGQAPRLLIYGASTRATGIPARFSGSGTEFTLTISSLQSEDSAVYYCQQYSNWPLTFGGGTKVEIK"]

# Insert a [MASK] token in the heavy chain (e.g., at the 10th position)
#masked_sequence = sequences_heavy[0][:9] + "[MASK]" + sequences_heavy[0][10:]
#paired_sequence = ' '.join(masked_sequence) + " [SEP] " + ' '.join(sequences_light[0])
paired_sequence = ' '.join(sequences_heavy[0]) + " [SEP] " + ' '.join(sequences_light[0])
# Tokenize the sequence
tokens = tokenizer(
    paired_sequence, 
    add_special_tokens=True, 
    padding="max_length", 
    return_tensors="pt"
)


# Select a position to mask (e.g., 10th token)
mask_position = 10  # Adjust as needed based on tokenized sequence
input_ids = tokens["input_ids"]

# Ensure the selected position is within bounds
if mask_position >= input_ids.shape[1]:
    raise ValueError("Selected mask position is out of bounds!")

# Replace the token at the chosen position with [MASK]
input_ids[0, mask_position] = tokenizer.mask_token_id

# Predict the masked token
with torch.no_grad():
    outputs = model(input_ids)
    predictions = outputs.logits

# Get the most probable token ID at the masked position
predicted_token_id = predictions[0, mask_position].argmax().item()
predicted_token = tokenizer.decode([predicted_token_id])

print(f"Predicted token for [MASK]: {predicted_token}")
