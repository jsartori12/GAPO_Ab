# GAPO-Ab: Genetic Algorithm for Protein Optimization - Antibody Version

GAPO-Ab is a specialized version of [GAPO](https://github.com/izzetbiophysicist/GAPO) specifically designed for antibody optimization. This tool implements a genetic algorithm approach to optimize antibody sequences for improved binding affinity and stability.

## Features

- Optimization of specific antibody regions:
  - CDRs (Complementarity-Determining Regions)
  - Framework regions
  - VH/VL interface
- Multiple mutation strategies using state-of-the-art language models:
  - ESM2 (Evolutionary Scale Modeling)
  - AbLang (antibody-specific language model)
  - Sapiens (human antibody repertoire-based model)
  - IgBert (paired antibody sequences model)
  - Random mutations
- Various objective functions for ΔGbind estimation:
  - Rosetta REF15 scoring function
  - Protein Binding Energy Estimator (PBEE)
  - Interface descriptors with Shape Complementarity

## Installation

```bash
# Clone the repository
git clone https://github.com/jsartori12/GAPO_Ab
cd GAPO_Ab

# Install dependencies (requirements.txt should be included in your repository)
pip install -r requirements.txt
```

## Dependencies

- PyRosetta
- Python 3.x
- Additional requirements will be listed in requirements.txt

## Algorithm Overview

1. **Initial Population Generation**
   - Random generation: Introduces random mutations at allowed positions
   - Optimized generation: Uses ESM2 model for naturally plausible sequences

2. **Structure Modeling**
   - Uses PyRosetta's Repack function for mutation modeling
   - Optimizes side chains within 8Å radius
   - Performs energy minimization using FastRelax

3. **Region Selection**
   - Implements ANARCI algorithm for CDR and framework identification
   - Interface detection between VH/VL domains (8Å distance criterion)

4. **Optimization Process**
   - Tournament selection
   - Crossover using binary mask
   - Multiple mutation strategies
   - Objective function evaluation

## Usage

[Add specific usage instructions and code examples here]

## Objective Functions

### 1. Rosetta ΔGbind
```
ΔGbind = ΔGcomplex - (ΔGAntibody + ΔGAntigen)
```

### 2. PBEE
Uses machine learning to predict binding free energy based on Rosetta structural descriptors.

### 3. Interface Analysis
```
ΔGInterface = Σ(InterfaceGiComplex) - (Σ(InterfaceGiAntibody) + Σ(InterfaceGiAntigen))
Final Score = ΔGInterface × Shape_Complementarity
```

## Citation

If you use GAPO-Ab in your research, please cite:
[Add citation information when available]

## License

[Add your chosen license]

## Contributors

[Add contributor information]

## Related Projects

- GAPO (Original version): [https://github.com/izzetbiophysicist/GAPO](https://github.com/izzetbiophysicist/GAPO)
- ESM2_NanoGEN: [https://github.com/jsartori12/ESM2_NanoGEN](https://github.com/jsartori12/ESM2_NanoGEN)

## Contact

[Add contact information]
