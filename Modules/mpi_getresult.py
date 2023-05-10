#!/usr/bin/python
# coding=utf-8
import os
import subprocess
import sys
import mpi4py.MPI as MPI


def get_filelist(dir):
    """
    获取目录中所有的文件
    :param dir:
    :return:
    """
    Filelist = []
    for home, dirs, files in os.walk(dir):
        for filename in files:
            # 文件名列表，包含完整路径
            if not filename == '.DS_Store':
                Filelist.append(os.path.join(home, filename))
    return Filelist


def get_pathlist(dir):
    """
    获取最底层的文件夹列表
    :param dir 根文件夹:
    :return 根文件夹下面所有最底层文件夹列表:
    """
    pathlist = []
    for home, dirs, files in os.walk(dir):
        if len(dirs) == 0 or len(files) != 0:
            for filename in files:
                if '.dok' in filename:
                    pathlist.append(home)
                    break
    return pathlist


def main(result_path):
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        pathlist = get_pathlist(result_path)

    else:
        pathlist = None

    pathlist = comm.bcast(pathlist, root=0)

    if rank == 0:
        data = []
        for i in range(1, size):
            tempData = comm.recv()
            data.extend(tempData)
        # 写入结果
        result_file_path = os.path.join(result_path, 'result.md')
        file_result = open(result_file_path, 'w')
        file_result.write(str(data))
        file_result.close()
    else:
        listTemp = []
        for fileName in get_filelist(pathlist[rank - 1]):
            try:
                ligandName = os.path.basename(fileName)
                receptorName = os.path.dirname(fileName).rpartition('/result/')[2]

                o = subprocess.Popen("sed -n 2p '" + fileName + "'|cut -d':' -f 3|sed 's/kcal\/mol//g'", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                oriscore = float(o.stdout.readlines()[0])
                n = subprocess.Popen("sed '/END/q' '" + fileName + "'|grep '^ATOM'|awk '{print $3}'|grep -v '^H'|wc -l", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                numatom = float(n.stdout.readlines()[0])
                lescore = oriscore / numatom

                _listTemp = [ligandName, receptorName, oriscore, numatom, lescore]
                listTemp.append(_listTemp)
            except BaseException as e:
                result_file_path = os.path.join(result_path, 'error.txt')
                file_result = open(result_file_path, 'a')
                file_result.write(str(e))
                file_result.write(str(fileName))
                file_result.write('\n')
                file_result.close()
                print(e)
                print(fileName)
                continue
        comm.send(listTemp, dest=0)


if __name__ == "__main__":
    try:
        result_path = sys.argv[1]
    except BaseException as e:
        # 测试
        result_path = '/Volumes/data/home/Developer/2021_06_29_00_13_37/result'
    finally:
        main(result_path)
