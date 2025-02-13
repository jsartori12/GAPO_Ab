#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 15:06:55 2025

@author: joao
"""

import sapiens
import argparse

parser = argparse.ArgumentParser(description="Process and complete protein sequences with masked tokens.")
parser.add_argument("--sequence", type=str, help="Protein sequence")
parser.add_argument("--H_L", type=str, help="Position to mask")
parser.add_argument("--jobid", type=str)
args = parser.parse_args()



to_model = args.sequence[0:110]

best = sapiens.predict_masked(to_model, args.H_L)

full_completed_sequence = best + args.sequence[len(to_model):]

with open("completed_sequence_"+args.jobid+".txt", "w") as f:
    f.write(full_completed_sequence)
    f.close()
