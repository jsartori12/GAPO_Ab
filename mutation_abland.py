#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 11:57:28 2025

@author: joao
"""

import ablang
import argparse

def complete_sequence(sequence, H_L):
    heavy_ablang = ablang.pretrained(H_L) # Use "light" if you are working with light chains
    heavy_ablang.freeze()


    seqs = []

    seqs.append(sequence)
    
    sequence = heavy_ablang(seqs, mode='restore')
    return sequence

parser = argparse.ArgumentParser(description="Process and complete protein sequences with masked tokens.")
parser.add_argument("--sequence", type=str, help="Protein sequence")
parser.add_argument("--H_L", type=str, help="Position to mask")
parser.add_argument("--jobid", type=str)
args = parser.parse_args()


to_model = args.sequence[0:110]


completed_sequence = complete_sequence(to_model, args.H_L)


full_completed_sequence = completed_sequence[0] + args.sequence[len(to_model):]


with open("completed_sequence_"+args.jobid+".txt", "w") as f:
    f.write(full_completed_sequence)
    f.close()
