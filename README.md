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

## Requirements

The program is originally written in Python 3.7
The following Python libraries are required to run it:
- rdkit
- numpy
- sklearn
- keras
- pandas
- bisect
- Levenshtein
- A backend to keras, such as theano, tensorflow or CNTK
