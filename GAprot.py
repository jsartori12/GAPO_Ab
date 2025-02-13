#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 03:52:23 2023

@author: lucas
"""

import os
from pyrosetta import Vector1
from genetic_algorithm_rosetta import genetic_algo, genetic_algo_sequence

from apt_function import *
import apt_function

from pyrosetta import *
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation

import numpy as np
from numpy.random import uniform
from random import sample
import random
import pandas as pd
import rosetta_descriptors


############ Running GAPO - Structure

# Initialize PyRosetta
pyrosetta.init()

# Create a scoring function using the "ref2015_cart.wts" weight set
scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

# Creates pose object from input PDB
starting_pose = pose_from_pdb('ab_trimed_relax.pdb')
# Relax the starting pose by packing and relaxing it iteratively for 3 times
scorefxn(starting_pose)

# Define a list of single-letter amino acid codes
gene_values = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
starting_pose_seq = [x for x in starting_pose.sequence()]
len(starting_pose_seq)


# Residues to mutate during optimization
CDRs = [412,
 413,
 414,
 415,
 416,
 431,
 432,
 433,
 434,
 435,
 436,
 437,
 438,
 439,
 440,
 441,
 442,
 443,
 444,
 445,
 446,
 479,
 480,
 481,
 482,
 483,
 484,
 485,
 486,
 487,
 488,
 489,
 194,
 195,
 196,
 197,
 198,
 199,
 200,
 201,
 202,
 203,
 204,
 220,
 221,
 222,
 223,
 224,
 225,
 226,
 259,
 260,
 261,
 262,
 263,
 264,
 265,
 266,
 267]
# Chain identifier for the fixed residues



# Generate an initial population of protein structures for optimization and return the list of residues to be locked during the evolution
init_population, list_fixed_index, listrosetta = apt_function.Generate_random_population(starting_pose = starting_pose, 
                                                             pop_size = 4,
                                                             fixed_residues_list = CDRs,
                                                             chains = ["C", "D"])

teste_pop = apt_function.Read_sequence("VHH_esm2_cdr_random_1.2.fasta")

teste_pop = teste_pop[:4]
new_init_pop = []

for i in teste_pop:
    temp = [*starting_pose_seq[0:171], *i]
    new_init_pop.append(temp)

Heavy_Light_chains = "C_D"

# Initiates GA object
GA = genetic_algo(pose=starting_pose, opt_direction='down',initial_population = new_init_pop, gene_values=gene_values, gene_type='discrete',
              vector_size=len(starting_pose_seq), pop_size=len(init_population), mutation_rate=0.9, segment_fluctuation=0,
              apt_function=apt_rosetta, selection_method='tournament', threads=False,
              convergence_threshold=0, n_cycles=3, tournament_cycles=int(np.round(len(init_population)/4)), tournament_size=4, benchmark=False,
              lista_fixed=list_fixed_index , crossing_over_type='mask', file_name="apt_rosetta_esm.txt", cpus  = 2, mutation_type = "ablang")

# Run the Genetic Algorithm
GA.execute()


    
# def unbind(pose, partners, scorefxn):
#     """
#     Simulate unbinding by applying a rigid body translation to a dummy pose and performing FastRelax.

#     Parameters:
#     - pose: PyRosetta Pose object representing the complex
#     - partners: Vector1 specifying the interacting partners
#     - scorefxn: Score function to evaluate the energy of the structure

#     Returns:
#     A tuple containing two Pose objects, one for the unbound state and one for the bound state.
#     """
#     #### Generates dummy pose to maintain original pose
    
#     pose_dummy = pose.clone()
#     pose_binded = pose.clone()
#     STEP_SIZE = 100
#     JUMP = 1
#     docking.setup_foldtree(pose_dummy, partners, Vector1([-1,-1,-1]))
#     trans_mover = rigid.RigidBodyTransMover(pose_dummy,JUMP)
#     trans_mover.step_size(STEP_SIZE)
#     trans_mover.apply(pose_dummy)
#     #pack_relax(pose_dummy, scorefxn)
#     #### Return a tuple containing:
#     #### Pose binded = [0] | Pose separated = [1]
#     return pose_binded , pose_dummy

# def Get_energy_contribution(pose, partner1, partner2, scorefxn):
    
#     partners = f"{partner1}_{partner2}"
    
#     binded, unbinded = unbind(pose, partners, scorefxn)
#     interface_residues = rosetta_descriptors.Get_interface_residues(pose = starting_posett, partner1 = partner1, partner2 = partner2)

#     scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
#     scorefxn.score(pose)
#     seq_size = len([x for x in pose.sequence()])
#     Residues = [residue.seqpos() for residue in pose]
#     df_byresidue_binded = pd.DataFrame(index=[1], columns=range(1, seq_size + 1))

#     residue_energies = [binded.energies().residue_total_energy(Residues[i-1]) for i in range(1, seq_size+1)]
#     df_byresidue_binded.iloc[0] = residue_energies

    
#     df_byresidue_unbinded = pd.DataFrame(index=[1], columns=range(1, seq_size + 1))

#     residue_energies = [unbinded.energies().residue_total_energy(Residues[i-1]) for i in range(1, seq_size+1)]
#     df_byresidue_unbinded.iloc[0] = residue_energies
    
    
    
#     return df_byresidue_binded, df_byresidue_unbinded

# teste1,teste2 = Get_energy_contribution(starting_posett, "A_D", scorefxn)

# binded, unbinded = unbind(starting_posett, "A_D", scorefxn)
# scorefxn(binded)
# scorefxn(unbinded)

# scorefxn(binded) - scorefxn(unbinded)



# binded, unbinded = unbind(starting_pose, "A_D", scorefxn)
# scorefxn(binded)
# scorefxn(unbinded)

# scorefxn(binded) - scorefxn(unbinded)



# def dG_interface(pose, to_unbind, partner1, partner2):
#     scorefxn(pose)

#     binded, unbinded = unbind(starting_pose, to_unbind, scorefxn)
#     interface_residues = rosetta_descriptors.Get_interface_residues(pose = pose, partner1 = partner1, partner2 = partner2)

#     scorefxn(binded)
#     scorefxn(unbinded)
#     # Calculate residue energies
#     seq_size = len(pose.sequence())
#     Residues = [residue.seqpos() for residue in pose]
    
#     # Binded residue energies
#     df_byresidue_binded = pd.DataFrame(index=[1], columns=range(1, seq_size + 1))
#     residue_energies_binded = [binded.energies().residue_total_energy(i) for i in Residues]
#     df_byresidue_binded.iloc[0] = residue_energies_binded
#     df_byresidue_binded.sum().sum()
#     # Unbinded residue energies
#     df_byresidue_unbinded = pd.DataFrame(index=[1], columns=range(1, seq_size + 1))
#     residue_energies_unbinded = [unbinded.energies().residue_total_energy(i) for i in Residues]
#     df_byresidue_unbinded.iloc[0] = residue_energies_unbinded
    
    
#     selected_columns = df_byresidue_binded.columns[interface_residues]
#     total_sum_binded = df_byresidue_binded[selected_columns].sum().sum()
    
    
#     selected_columns_unb = df_byresidue_unbinded.columns[interface_residues]
#     total_sum_unbinded = df_byresidue_unbinded[selected_columns_unb].sum().sum()
    
#     ddG_interface = total_sum_binded - total_sum_unbinded
    
#     testeposesss, testedecriptorssss = Get_descriptors(pose = starting_pose, partner1 = "CD", partner2 = "A")
    
#     sc_value = testedecriptorssss["ifa_sc_value"][0]
    
#     apt_dg_sc = ddG_interface * sc_value
    
#     return apt_dg_sc


# testedg_interface = dG_interface(pose = starting_pose, to_unbind = "A_D", partner1 = "CD", partner2 = "A")


# binded.dump_pdb("bindedteste.pdb")
# unbinded.dump_pdb("unbindedteste.pdb")

# residue_index = 10  # Example residue
# print("Binded residue energy:", binded.energies().residue_total_energy(residue_index))
# print("Unbinded residue energy:", unbinded.energies().residue_total_energy(residue_index))

# teste

# energies = unbinded.energies()


# podebinded = pose_from_pdb('bindedteste.pdb')
# podeunbinded= pose_from_pdb('unbindedteste.pdb')

# scorefxn(podebinded)
# scorefxn(podeunbinded)


# def Energy_contribution(pose, by_term):
#     """
#     Calculate and analyze the energy contributions of different terms for each residue in a given protein pose.

#     Parameters:
#     - pose: PyRosetta Pose object representing the protein structure.
#     - by_term: Boolean flag indicating whether to analyze energy contributions by term (True) or by residue (False).

#     Returns:
#     - DataFrame: If by_term is True, returns a DataFrame containing energy contributions by term for each residue.
#                  If by_term is False, returns a DataFrame containing energy contributions by residue.
#     """

#     # List of energy terms from ref2015_cart
#     listadosdicts = [fa_atr, fa_rep, fa_sol, fa_intra_rep, fa_intra_sol_xover4,
#                 lk_ball_wtd, fa_elec, hbond_sr_bb, hbond_lr_bb, hbond_bb_sc, hbond_sc, dslf_fa13,
#                 omega, fa_dun, p_aa_pp, yhh_planarity, ref, rama_prepro, cart_bonded]
    
#     # Create a score function using ref2015_cart weights
#     scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
#     weights = scorefxn.weights()
    
#     # Set up energy method options to decompose hbond terms
#     emopts = EnergyMethodOptions(scorefxn.energy_method_options())
#     emopts.hbond_options().decompose_bb_hb_into_pair_energies(True)
#     scorefxn.set_energy_method_options(emopts)
    
#     # Calculate energy scores for the given pose
#     scorefxn.score(pose)
    
#     # Check if the user wants to analyze energy contributions by term
#     if by_term == True:
#         # Initialize a dictionary for storing data
#         dasd = {'Protein': [], 'Sequence': []}
        
#         # Get all residues' pose index from pose
#         Residues = [residue.seqpos() for residue in pose]
#         dasd['Protein'].append("WT")
        
#         # Populate the dictionary with energy contributions for each residue and term
#         for posi in Residues:
#             for i in range(len(listadosdicts)):
#                 term_key = '{}-%s'.format(posi) % listadosdicts[i]
#                 dasd[term_key] = []
#                 dasd[term_key].append(pose.energies().residue_total_energies(posi)[listadosdicts[i]])
#         dasd['Sequence'].append(pose.sequence())
        
#         # Create a DataFrame from the dictionary
#         df2 = pd.DataFrame(dasd)
        
#         # Create a DataFrame with energy terms and their respective weights
#         weights_by_term = pd.DataFrame(index=range(1, len(listadosdicts)+1), columns=range(0, 2))
#         weights_by_term.iloc[:, 0] = listadosdicts
#         list_weights = [1, 0.55, 1, 0.005, 1, 1, 1, 1, 1, 1, 1, 1.25, 0.4, 0.7, 0.6, 0.625, 1, 0.45, 0.5]
#         weights_by_term.iloc[:, 1] = list_weights
        
#         # Apply weights to each term in the energy contribution DataFrame
#         for i in range(len(weights_by_term)):
#             list_to_change = df2.filter(like=str(weights_by_term.iloc[i, 0])).columns
#             df2[list_to_change] = df2[list_to_change] * weights_by_term.iloc[i, 1]
        
#         return df2
#     else:
#         # If not analyzing by term, create a DataFrame with energy contributions by residue
#         seq_size = len([x for x in pose.sequence()])
#         Residues = [residue.seqpos() for residue in pose]
#         df_byresidue = pd.DataFrame(index=range(1, 2), columns=range(1, seq_size+1))
        
#         for i in range(1, len(df_byresidue.columns)+1):
#             df_byresidue.iloc[0, i-1] = pose.energies().residue_total_energy(Residues[i-1])
        
#         return df_byresidue
# #imports from pyrosetta
# from mimetypes import init
# from pyrosetta import *
# from pyrosetta.teaching import *
# from rosetta.core.kinematics import MoveMap
# from rosetta.core.kinematics import FoldTree
# from rosetta.core.pack.task import TaskFactory
# from rosetta.core.pack.task import operation
# from rosetta.core.simple_metrics import metrics
# from rosetta.core.select import residue_selector as selections
# from rosetta.core import select
# from rosetta.core.select.movemap import *
# from rosetta.protocols import minimization_packing as pack_min
# from rosetta.protocols import relax as rel
# from rosetta.protocols.antibody.residue_selector import CDRResidueSelector
# from rosetta.protocols.antibody import *
# from rosetta.protocols.loops import *
# from rosetta.protocols.relax import FastRelax
# from pyrosetta.rosetta.protocols.docking import setup_foldtree
# from pyrosetta.rosetta.protocols import *
# from rosetta.core.scoring.methods import EnergyMethodOptions
# all_binded = Energy_contribution(binded, True)
# all_binded_values = all_binded.iloc[:, 2:]
# all_unbinded = Energy_contribution(unbinded, True)
# all_unbinded_values = all_unbinded.iloc[:, 2:]

# sum_columns_binded = all_binded.iloc[:, 2:].sum(axis=1)
# sum_columns_unbinded = all_unbinded.iloc[:, 2:].sum(axis=1)


# teste

# column_indices = []
# for residue in teste:
#     # Calculate the starting column index for the residue
#     start_col = 1 + (residue - 1) * 19  # Each residue has 12 energy terms
#     # Add all 12 columns for the residue
#     column_indices.extend(range(start_col, start_col + 12))

# # Extract the relevant columns
# extracted_columns = all_binded_values.iloc[:, column_indices]
# extracted_columnsunbid = all_unbinded_values.iloc[:, column_indices]

# extracted_columns.sum().sum()
# extracted_columnsunbid.sum().sum()





















