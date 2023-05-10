from mol_utils import drop_salt
from rdkit import Chem

# Read a file containing SMILES
# The file should be a .smi or a .csv where the first column should contain a SMILES string
def read_file(file_name, drop_first=True):
    
    molObjects = []

    with open(file_name) as f:
        for l in f:
            if drop_first:
                drop_first = False
                continue

            l = l.strip().split(",")[0]  #strip()表示删除掉数据中的换行符，split（‘，’）则是数据中遇到‘,’ 就隔开。把每行的每个字符一个个分开，变成一个list
            smi = drop_salt(l.strip())
            molObjects.append(Chem.MolFromSmiles(smi))

    return molObjects

