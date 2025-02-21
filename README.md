GAPO - Genetic Algorithm for Protein Optimization

Overview

GAPO (Genetic Algorithm for Protein Optimization) is an in silico genetic algorithm designed to optimize protein sequences based on structural and energetic evaluations. It allows for directed evolution simulations by iterating through mutation, recombination, and selection processes to optimize stability, binding affinity, or other protein properties.

The framework is implemented using PyRosetta for structural modeling and includes multiple objective functions based on both physics-based and machine learning models.

Variants:

GAPO: A general-purpose version for protein optimization, allowing stability and sequence probability optimization using ESM2.

GAPO-ab: A specialized version for antibody optimization, incorporating modules for CDR and framework detection, as well as VH/VL interface optimization. It also integrates antibody-specific language models for guided mutations.

Algorithm Workflow

Population Generation

Random Generation: Sequences are initialized with random mutations in permitted positions, allowing for a broad search across sequence space.

Optimized Generation: Uses ESM2_NanoGEN, a model leveraging Evolutionary Scale Modeling (ESM), to introduce mutations informed by evolutionary constraints, ensuring biologically plausible sequences.

Structural Modeling

Each sequence is modeled using PyRosetta's Repack function, mutating residues and optimizing side chains within an 8 Å radius.

Energy minimization is performed with FastRelax to refine structures before evaluation.

Objective Functions
The algorithm evaluates candidates using different scoring functions:

Rosetta ∆Gbind: Uses the REF15 score function to estimate binding energy as:


Protein Binding Energy Estimator (PBEE): Uses a machine learning model trained on Rosetta-derived structural descriptors to predict binding energy (PBEE GitHub).

Rosetta Interface Descriptors: Uses InterfaceAnalyzerMover to extract interface-specific energy contributions and shape complementarity (SC) scores:

The objective function is defined as the product , emphasizing both binding energy and shape complementarity.

Optimization Process

Tournament Selection: Four candidates are randomly selected, and the two with the best binding energy (∆Gbind) proceed.

Crossover: Two-parent recombination uses a binary mask to generate offspring by mixing sequence positions.

Mutation Strategies:

Random: Unbiased random mutations for broad exploration.

ESM2: Context-aware mutations predicted by ESM2.

AbLang: Mutations informed by an antibody-specific language model trained on unpaired sequences.

Sapiens: Uses human antibody repertoires for humanization-driven mutations.

IgBert: Leverages paired antibody chains for optimized VH/VL interactions.

Iterative Optimization

The process repeats through multiple generations until convergence or a predefined stopping criterion is met.

GAPO vs. GAPO-ab

Feature

GAPO

GAPO-ab

Stability Optimization (∆Gfold)

✓

✓

ESM2 Sequence Probability Optimization

✓

✓

CDR/Framework Identification



✓

VH/VL Interface Optimization



✓

Antibody-Specific Mutation Models



✓

Affinity Optimization (∆Gbind)

✓

✓

ML-based Binding Prediction (PBEE)

✓

✓

Rosetta Interface Score Optimization

✓

✓

Installation & Usage

Dependencies

Python 3.8+

PyRosetta

ESM2

ANARCI (for antibody numbering)

PBEE (optional for ML-based scoring)

Running GAPO

python run_gapo.py --config config.json

Running GAPO-ab

python run_gapo_ab.py --config config_ab.json

References

Lin et al. 2023, Evolutionary Scale Modeling (ESM)

Hie et al. 2024, Antibody Optimization Using ML

Chaudhury, Lyskov, and Gray 2010, PyRosetta

For more details, visit our repositories:

GAPO: GitHub

GAPO-ab: GitHub

