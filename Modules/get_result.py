#!/usr/bin/python
# coding=utf-8
import datetime
import os
import sys
import logging

logging.basicConfig(
    filename='get_result.log',
    level=logging.DEBUG,
    format=' %(asctime)s - %(levelname)s - [文件：%(filename)s 第%(lineno)d行 方法：%(funcName)s] - %(message)s'
)
logging.captureWarnings(True)


def get_filelist(dir):
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


def main(currentPath):
    logging.info('=' * 100)
    logging.info(str((datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))))
    logging.info('=' * 100)

    logging.info(currentPath)

    pahtList = get_pathlist(currentPath)

    pathSize = len(pahtList)
    if pathSize == 0:
        logging.error('没有生成结果')
        return '没有生成结果！'
    pathSize = pathSize + 1

    # 高研院
    ppn = 36

    # 即墨超算
    # ppn = 28
    nodes = int(pathSize / ppn) + 1

    qsubFileName = os.path.join(currentPath, 'get_result.pbs')
    file_qsub = open(qsubFileName, 'w')
    file_qsub.write('#pbs.sh' + '\n')
    file_qsub.write('#!/bin/bash' + '\n')
    file_qsub.write('#PBS -N ' + 'GetResult\n')
    file_qsub.write('#PBS -j oe' + '\n')
    file_qsub.write('#PBS -l nodes=' + str(nodes) + ':ppn=' + str(ppn) + '\n')
    file_qsub.write('#PBS -l walltime=999:00:00' + '\n')
    file_qsub.write('#PBS -q com' + '\n')
    # file_qsub.write('#PBS -q q_csywz' + '\n')
    file_qsub.write('#PBS -V' + '\n')
    file_qsub.write('#PBS -S /bin/bash' + '\n')
    # file_qsub.write('#PBS -o ' + currentPath + '/logs/' + qsubName + '.log' + '\n')
    file_qsub.write('\n')
    file_qsub.write('date' + '\n')
    file_qsub.write('mpiexec -n ' + str(pathSize) + ' python ' + os.path.join(currentPath.rpartition('/result')[0], 'mpi_getresult.py ') + currentPath + '\n')
    file_qsub.write('date' + '\n')
    file_qsub.close()

    logging.debug('生成脚本文件：' + qsubFileName)
    # 更改脚本执行目录
    os.chdir(currentPath)
    # 提交获取结果任务


    # TODO 高研院
    qsub_command = '/usr/bin/qsub ' + qsubFileName

    # TODO 即墨超算
    # qsub_command = '/home/SystemSoftware/tsce/torque6/bin/qsub ' + qsubFileName

    os.system(qsub_command)

if __name__ == "__main__":
    try:
        currentPath = sys.argv[1]
    except BaseException as e:
        # 测试
        currentPath = '/Volumes/data/home/Downloads/result'
    finally:
        main(currentPath)
