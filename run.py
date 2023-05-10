import sys
sys.path.insert(0, './Modules/')

import numpy as np

from Modules.file_reader import read_file
from Modules.mol_utils import get_fragments
from Modules.build_encoding import get_encodings, encode_molecule, decode_molecule, encode_list, save_decodings
from Modules.models_alter_transformer import build_models
from Modules.training import train
from Modules.rewards import clean_good
from rdkit import rdBase
import logging
logging.getLogger().setLevel(logging.INFO)
rdBase.DisableLog('rdApp.error')


def main(fragment_file, lead_file):
    # 将smi形式转换为mol形式的列表
    fragment_mols = read_file(fragment_file)
    print(fragment_mols[0])
    lead_mols = read_file(lead_file)
    # fragment将所有分子存入
    fragment_mols += lead_mols
    #print(fragment_mols)   rdkit.Chem.rdchem.Mol object at 0x0000013B02F4DE70
    print("Read {} molecules for fragmentation library", len(fragment_mols))   #logging.info（）  输出日志的信息
    print("Read {} lead moleculs", len(lead_mols))

    fragments, used_mols = get_fragments(fragment_mols)
    print("Num fragments: {}".format(len(fragments)) )
    print("Total molecules used: {}", len(used_mols))
    assert len(fragments)
    assert len(used_mols)
    encodings, decodings = get_encodings(fragments)
    save_decodings(decodings)
    # #print(decodings)   '0000001': <rdkit.Chem.rdchem.Mol object at 0x0000013B02F4DE70>
    logging.info("Saved decodings")
    
    lead_mols = np.asarray(fragment_mols[-len(lead_mols):])[used_mols[-len(lead_mols):]]
    # #print(lead_mols) <rdkit.Chem.rdchem.Mol object at 0x0000022C01F30490>   47个
    print(len(lead_mols))
    X = encode_list(lead_mols, encodings)
    print(X)
    print(X.shape)
    # logging.info("Building models")
    # #print(X.shape[0]) 47         print(X.shape[1])  最大12个片段     print(X.shape[2])  8 （一个片段为8位）
    print(X.shape[1:])
    actor, critic = build_models(X.shape[1:])
    # print(actor.summary())

    # #print(X)

    # logging.info("Training")
    history = train(X, actor, critic, decodings)

    # logging.info("Saving")
    np.save("History/history.npy", history)




if __name__ == "__main__":

    fragment_file = "Data/molecules.smi"
    lead_file = "Data/dopamineD4props.csv"


    if len(sys.argv) > 1:
        fragment_file = sys.argv[1]

    if len(sys.argv) > 2:
        lead_file = sys.argv[2]

    main(fragment_file, lead_file)


