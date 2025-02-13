#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 11:05:09 2024

@author: joao
"""

from pyrosetta import rosetta
from pyrosetta.rosetta.protocols import *
from rosetta.core.pack.task import TaskFactory, operation
from rosetta.core.kinematics import MoveMap, FoldTree
from rosetta.core.simple_metrics import metrics
from rosetta.core.select import residue_selector as selections
from rosetta.core import select
from rosetta.core.select.movemap import *
from rosetta.protocols import minimization_packing as pack_min
from rosetta.protocols import relax as rel
from rosetta.protocols.antibody.residue_selector import CDRResidueSelector
from rosetta.protocols.antibody import *
from rosetta.protocols.loops import *
from rosetta.protocols.relax import FastRelax
from rosetta.protocols.docking import setup_foldtree
from pyrosetta import Vector1
from rosetta.core.scoring.methods import EnergyMethodOptions

import pyrosetta
import pandas as pd
import os
import sys

def jd2_format(pdbfile, basename, outdir):
    pyrosetta.init(extra_options="-corrections::beta_nov16 true -ignore_unrecognized_res -output_pose_energies_table false -renumber_pdb")
    pose = rosetta.core.import_pose.pose_from_file(pdbfile)
    pose.dump_pdb(f'{outdir}/{basename}_jd2_0001.pdb')

# Define a função de minimização
def minimize(pose, scorefxn, minimizer_type):
    # Cria um MoveMap dependendo do tipo de minimização
    movemap = pyrosetta.rosetta.core.kinematics.MoveMap()
    if minimizer_type == 'minmover1':
        movemap.set_bb(False)
        movemap.set_chi(True)
        tolerance = 0.0001
    elif minimizer_type == 'minmover2':
        movemap.set_bb(True)
        movemap.set_chi(True)
        tolerance = 0.0001

    # Cria o MinMover
    min_mover = rosetta.protocols.minimization_packing.MinMover()
    min_mover.score_function(scorefxn)
    min_mover.max_iter(50000)
    min_mover.tolerance(tolerance)
    min_mover.cartesian(False)
    min_mover.movemap(movemap)

    # Aplica a minimização à pose
    min_mover.apply(pose)
    
def unbind(pose, partners, scorefxn):
    """
    Simulate unbinding by applying a rigid body translation to a dummy pose and performing FastRelax.

    Parameters:
    - pose: PyRosetta Pose object representing the complex
    - partners: Vector1 specifying the interacting partners
    - scorefxn: Score function to evaluate the energy of the structure

    Returns:
    A tuple containing two Pose objects, one for the unbound state and one for the bound state.
    """
    #### Generates dummy pose to maintain original pose
    
    pose_dummy = pose.clone()
    pose_binded = pose.clone()
    STEP_SIZE = 100
    JUMP = 1
    docking.setup_foldtree(pose_dummy, partners, Vector1([-1,-1,-1]))
    trans_mover = rigid.RigidBodyTransMover(pose_dummy,JUMP)
    trans_mover.step_size(STEP_SIZE)
    trans_mover.apply(pose_dummy)
    #pack_relax(pose_dummy, scorefxn)
    #### Return a tuple containing:
    #### Pose binded = [0] | Pose separated = [1]
    return pose_binded , pose_dummy

def Get_energy_per_term(pose, scorefxn):
    # Get and display the individual weighted energy terms
    score_types = rosetta.core.scoring.ScoreType
    weights = scorefxn.weights()
    energy_map = pose.energies().total_energies()  # Get total energies for the pose
    energy_terms = {}
    
    # Loop through all score types
    for score_type in rosetta.core.scoring.ScoreType.__members__.values():
        weight = weights[score_type]
        if weight != 0:  # Only include terms with non-zero weights
            term_value = energy_map[score_type] * weight
            energy_terms[score_type.name] = term_value
    return energy_terms

def Interaction_energy_metric(pose, scorefxn, partner1, partner2):
    # Create the OrResidueSelector for partner1
    partner1_selector = rosetta.core.select.residue_selector.OrResidueSelector()
    
    # Add ChainSelectors for each chain in partner1
    for chain in partner1:
        chain_selector = rosetta.core.select.residue_selector.ChainSelector(chain)
        partner1_selector.add_residue_selector(chain_selector)
    
    # Create the OrResidueSelector for partner2
    partner2_selector = rosetta.core.select.residue_selector.OrResidueSelector()
    
    # Add ChainSelectors for each chain in partner2
    for chain in partner2:
        chain_selector = rosetta.core.select.residue_selector.ChainSelector(chain)
        partner2_selector.add_residue_selector(chain_selector)
    
    # Initialize InteractionEnergyMetric
    interaction_energy_metric = rosetta.core.simple_metrics.metrics.InteractionEnergyMetric()
    interaction_energy_metric.set_scorefunction(scorefxn)
    
    # Set both residue selectors
    interaction_energy_metric.set_residue_selectors(partner1_selector, partner2_selector)
    
    # Calculate the interaction energy metric
    interaction_energy = interaction_energy_metric.calculate(pose)
    
    # Return the result
    return interaction_energy

def Contact_molecular_surface(pose, partner1, partner2):
    # Create the OrResidueSelector for partner1
    partner1_selector = rosetta.core.select.residue_selector.OrResidueSelector()
    
    # Add ChainSelectors for each chain in partner1
    for chain in partner1:
        chain_selector = rosetta.core.select.residue_selector.ChainSelector(chain)
        partner1_selector.add_residue_selector(chain_selector)
    
    # Create the OrResidueSelector for partner2
    partner2_selector = rosetta.core.select.residue_selector.OrResidueSelector()
    
    # Add ChainSelectors for each chain in partner2
    for chain in partner2:
        chain_selector = rosetta.core.select.residue_selector.ChainSelector(chain)
        partner2_selector.add_residue_selector(chain_selector)
    
    cms_filter = rosetta.protocols.simple_filters.ContactMolecularSurfaceFilter()
    cms_filter.selector1(partner1_selector)
    cms_filter.selector2(partner2_selector)
    cms_filter.distance_weight(0.5)
    cms_filter.set_user_defined_name("cms")
    cms_filter.apply(pose)
    return cms_filter.score(pose)

def Interface_analyzer_mover(pose, partner1, partner2):
    partners = f"{partner1}_{partner2}"
    ifa_mover = rosetta.protocols.analysis.InterfaceAnalyzerMover(partners)
    ifa_mover.set_use_tracer(True)
    ifa_mover.set_compute_packstat(True)
    ifa_mover.set_scorefile_reporting_prefix("ifa")
    ifa_mover.apply(pose)
    interface_data = ifa_mover.get_all_data()
    ifa_mover.add_score_info_to_pose(pose)
    score_map = pose.scores
    # Filter the elements where the key starts with "ifa"
    ifa_data = {score_name: value for score_name, value in score_map.items() if score_name.startswith("ifa")}
    return ifa_data

def Get_interface_selector(pose, partner1, partner2):
    partners  = f"{partner1}_{partner2}"
    ifa_mover = rosetta.protocols.analysis.InterfaceAnalyzerMover(partners)
    ifa_mover.set_use_tracer(True)
    ifa_mover.set_compute_packstat(True)
    ifa_mover.set_scorefile_reporting_prefix("ifa")
    ifa_mover.apply(pose)
    interface_data = ifa_mover.get_all_data()
    residues_in_interface = []
    for i in range(1, len(interface_data.interface_residues[1]) + 1):
        if interface_data.interface_residues[1][i]:
            residues_in_interface.append(i)
    residue_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    residue_selector.set_index(','.join(map(str, residues_in_interface)))
    return residue_selector

def Get_interface_residues(pose, partner1, partner2):
    partners  = f"{partner1}_{partner2}"
    ifa_mover = rosetta.protocols.analysis.InterfaceAnalyzerMover(partners)
    ifa_mover.set_use_tracer(True)
    ifa_mover.set_compute_packstat(True)
    ifa_mover.set_scorefile_reporting_prefix("ifa")
    ifa_mover.apply(pose)
    interface_data = ifa_mover.get_all_data()
    residues_in_interface = []
    for i in range(1, len(interface_data.interface_residues[1]) + 1):
        if interface_data.interface_residues[1][i]:
            residues_in_interface.append(i)
    return residues_in_interface

def Get_energy_contribution(pose, partners, scorefxn):
    
    binded, unbinded = unbind(pose, partners, scorefxn)
    
    scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
    scorefxn.score(pose)
    seq_size = len([x for x in pose.sequence()])
    Residues = [residue.seqpos() for residue in pose]
    df_byresidue_binded = pd.DataFrame(index=range(1, 2), columns=range(1, seq_size+1))
    
    for i in range(1, len(df_byresidue_binded.columns)+1):
        df_byresidue_binded.iloc[0, i-1] = binded.energies().residue_total_energy(Residues[i-1])
    
    df_byresidue_unbinded = pd.DataFrame(index=range(1, 2), columns=range(1, seq_size+1))
    
    for i in range(1, len(df_byresidue_unbinded.columns)+1):
        df_byresidue_unbinded.iloc[0, i-1] = unbinded.energies().residue_total_energy(Residues[i-1])
    
    
    
    return df_byresidue_binded, df_byresidue_unbinded



def Get_descriptors(pose, partner1, partner2):

    pyrosetta.init(extra_options="\
    -mute core \
    -mute basic \
    -ex1 -ex2 -ex1aro -ex2aro \
    -use_input_sc \
    -flip_HNQ \
    -no_optH false \
    -corrections::beta_nov16 true \
    -output_pose_energies_table false")
    #pose = rosetta.core.import_pose.pose_from_file(pdb)
    pose = pose
    
    # Create the beta_nov16 score function
    scorefxn = rosetta.core.scoring.ScoreFunctionFactory.create_score_function("beta_nov16")
    scorefxn(pose)
    # # Set up energy method options to decompose hbond terms
    # emopts = EnergyMethodOptions(scorefxn.energy_method_options())
    # emopts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    # scorefxn.set_energy_method_options(emopts)
     
 # Calculate energy scores for the given pose
    scorefxn.score(pose)
    # minimize(pose=pose, scorefxn=scorefxn, minimizer_type="minmover1")
    # minimize(pose=pose, scorefxn=scorefxn, minimizer_type="minmover2")
    # scorefxn(pose)

    #per_term  = Get_energy_per_term(pose, scorefxn)
    #testeIE   = Interaction_energy_metric(pose = pose, scorefxn = scorefxn, partner1 = partner1, partner2 = partner2)
    #testeCMS  = Contact_molecular_surface(pose = pose, partner1 = partner1, partner2 = partner2)
    all_terms = Interface_analyzer_mover(pose  = pose, partner1 = partner1, partner2 = partner2)
    #all_terms = per_term | testeIFA

    #all_terms["cms"] = testeCMS
    #all_terms["interaction_energy"] = testeIE
    #all_terms["total_score"] = scorefxn(pose)
    all_terms_df = pd.DataFrame([all_terms])
    all_terms_df.rename(columns={
        "ifa_dG_separated/dSASAx100":"ifa_dG_separated_dSASAx100", 
        "ifa_dG_cross/dSASAx100":"ifa_dG_cross_dSASAx100"}, inplace=True)
    #new_apt = all_terms_df["ifa_dG_separated"].value / all_terms_df["ifa_sc_value"].value
    #print(f"New apt_function: {new_apt}")
    return pose, all_terms_df

def dG_interface(pose, to_unbind, partner1, partner2):
    
    scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
    scorefxn(pose)

    binded, unbinded = unbind(pose, to_unbind, scorefxn)
    interface_residues = Get_interface_residues(pose = pose, partner1 = partner1, partner2 = partner2)

    scorefxn(binded)
    scorefxn(unbinded)
    
        
    # Set up energy method options to decompose hbond terms
    emopts = EnergyMethodOptions(scorefxn.energy_method_options())
    emopts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    scorefxn.set_energy_method_options(emopts)
    
    # Calculate energy scores for the given pose
    scorefxn.score(pose)
    
    # Calculate residue energies
    seq_size = len(pose.sequence())
    Residues = [residue.seqpos() for residue in pose]
    weights =  local_scorefxn.weights() 

    emopts = core.scoring.methods.EnergyMethodOptions( local_scorefxn.energy_method_options() ) 
    emopts.hbond_options().decompose_bb_hb_into_pair_energies( True )
    local_scorefxn.set_energy_method_options( emopts )
    # Binded residue energies
    df_byresidue_binded = pd.DataFrame(index=[1], columns=range(1, seq_size + 1))
    residue_energies_binded = [binded.energies().residue_total_energy(i) for i in Residues]
    df_byresidue_binded.iloc[0] = residue_energies_binded
    df_byresidue_binded.sum().sum()
    # Unbinded residue energies
    df_byresidue_unbinded = pd.DataFrame(index=[1], columns=range(1, seq_size + 1))
    residue_energies_unbinded = [unbinded.energies().residue_total_energy(i) for i in Residues]
    df_byresidue_unbinded.iloc[0] = residue_energies_unbinded
    
    
    selected_columns = df_byresidue_binded.columns[interface_residues]
    total_sum_binded = df_byresidue_binded[selected_columns].sum().sum()
    
    
    selected_columns_unb = df_byresidue_unbinded.columns[interface_residues]
    total_sum_unbinded = df_byresidue_unbinded[selected_columns_unb].sum().sum()
    
    ddG_interface = total_sum_binded - total_sum_unbinded
    
    testeposesss, testedecriptorssss = Get_descriptors(pose = pose, partner1 = "CD", partner2 = "A")
    
    sc_value = testedecriptorssss["ifa_sc_value"][0]
    
    apt_dg_sc = ddG_interface * sc_value
    
    return apt_dg_sc
