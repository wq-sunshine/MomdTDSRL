# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 16:00:57 2021

@author: b519
"""


import os
path ='D:/generatedMolecules/DeepFMPO/deep/ledock/Mpro/Mpro'
def get_filelist(dir):
    Filelist = []
    for home, dirs, files in os.walk(path):
        for filename in files:
            # 文件名列表，包含完整路径
            Filelist.append(os.path.join(home, filename))
            # # 文件名列表，只包含文件名
            # Filelist.append( filename)
            
    return Filelist
if __name__ =="__main__":
    Filelist = get_filelist(dir)
    print(len( Filelist))
    #old_str = "pro.pdb"
    #old_str = "ligandsn"
    for file in  Filelist :
        if file.endswith('_ledock.in'):
            #file_data = ""
            p1,f1 = os.path.split(file)
            p2,f2 = os.path.split(p1)
            #new_str = os.path.join(p1,f2+"_pro.pdb")
            #new_str = "/home/b519/lyy/lyy2/DeepFMPO/DeepFMPO-master/3.6/ledock/ligand0"
            #with open(file, "r", encoding="utf-8") as f:
                #for line in f:
                    #if old_str == line:
                        #line = line.replace(old_str,new_str)
                    #file_data += line
            #with open(file,"w",encoding="utf-8") as f:
                #f.write(file_data)
            #os.system('cd '+p1)
            #os.system('/home/b519/lyy/lyy2/DeepFMPO/DeepFMPO-master/3.6/ledock/ledock_go '+file)
            os.system('sh dock.sh D:/generatedMolecules/DeepFMPO/deep/ledock/result '+file+' '+f2)
            print(file)
        
