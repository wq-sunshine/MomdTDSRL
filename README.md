# MomdTDSRL
Code accompanying the paper "Multi-objective molecular design based on transformer encoder, docking simulations and reinforcement learning".

## Instructions

To run the main program on the same data as described in the paper just run:
```sh
python run.py
```
It is also possible to run the program on a custom set of lead molecules and/or fragments. 
```sh
python run.py fragment_molecules.smi lead_file.smi
```
Molecules that are genertated during the process can be vied by running:
```sh
python Show_Epoch.py n
```
where `n` is the epoch that should be viewed. This shows two columns of molecules. The first column contains the original lead molecule, while the second column contains modified molecules.
Any global parameters can be changed by changing them in the file "Modules/global_parameters.py"

## Data collection
The virtual fragment combinatorial library: fragmenting compounds with 9823 Mpro inhibitory activity from the ChEMBL database.
The initial lead compound molecular group: 175 lead compounds.
In the pretraining process of the perceptron classifier, the positive and negative sample data were from the dekois 2.0 benchmark data set library provided by Tubingen University that contained the SARS coronavirus 3CL protease data set.

## Requirements

The program is originally written in Python 3.7
The following Python libraries are required to run it:
- rdkit
- ledock 1.0
- numpy 1.19.5
- sklearn
- matplotlib 3.1.0
- pandas 0.24.2
- scikit-learn 0.24.1
- Levenshtein 0.16.0
- A backend to keras 2.4.3, such as tensorflow 2.2.2
