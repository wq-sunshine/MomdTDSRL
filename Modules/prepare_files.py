#!/usr/bin/python
# coding=utf-8
import datetime
import os
import shutil
import time

import logging
import zipfile

logging.basicConfig(
    filename='prepare_files.log',
    level=logging.DEBUG,
    format=' %(asctime)s - %(levelname)s - [文件：%(filename)s 第%(lineno)d行 方法：%(funcName)s] - %(message)s'
)
logging.captureWarnings(True)

'''
生产脚本文件
'''


def receptor_qsub(nodes, ppn, count, currentPath, ledockConfigPathList):
    start_time = time.time()
    qsubPath = os.path.join(currentPath, 'qsub_pbs')
    qsubPathList = []
    ledockNameList = []
    for ledockConfigPath in ledockConfigPathList:
        ledockName = os.path.basename(ledockConfigPath)
        ledockNameList.append(ledockName)
    for ledcokNameSplitList in [ledockNameList[i:i + count] for i in range(0, len(ledockNameList), count)]:
        qsubFileName = os.path.join(qsubPath, '^'.join(ledcokNameSplitList) + '.pbs')
        file_qsub = open(qsubFileName, 'w')
        file_qsub.write('#pbs.sh' + '\n')
        file_qsub.write('#!/bin/bash' + '\n')
        file_qsub.write('#PBS -N ' + '^'.join(ledcokNameSplitList) + '\n')
        file_qsub.write('#PBS -j oe' + '\n')
        file_qsub.write('#PBS -l nodes=' + str(nodes * len(ledcokNameSplitList)) + ':ppn=' + str(ppn) + '\n')
        file_qsub.write('#PBS -l walltime=999:00:00' + '\n')
        file_qsub.write('#PBS -q com' + '\n')
        # file_qsub.write('#PBS -q q_csywz' + '\n')
        file_qsub.write('#PBS -V' + '\n')
        file_qsub.write('#PBS -S /bin/bash' + '\n')
        # file_qsub.write('#PBS -o ' + currentPath + '/logs/' + qsubName + '.log' + '\n')
        file_qsub.write('\n')
        file_qsub.write('date' + '\n')
        for ledockName in ledcokNameSplitList:
            ledockpathlist = list(filter(lambda f: str(f).endswith(ledockName), ledockConfigPathList))
            for ledockpath in ledockpathlist:
                file_qsub.write('mpiexec -n ' + str(nodes * ppn) + ' python ' + os.path.join(currentPath, 'mpi_ledock.py ')
                                + currentPath + ' ' + os.path.join(currentPath, 'result ') + ledockpath + '\n')
        file_qsub.write('date' + '\n')
        file_qsub.close()
        qsubPathList.append(qsubFileName)

    used_time = time.time() - start_time
    logging.info("生成脚本文件所使用时间：%s" % (str(datetime.timedelta(seconds=used_time))))
    return qsubPathList


'''
调整配置文件
'''


def Split_file_pdb(currentPath, receptorpath, ligandFileListNew):
    start_time = time.time()
    ligandNum = len(ligandFileListNew)
    # ledock 配置文件目录
    ledockinPath = os.path.join(currentPath, 'ledockin')
    if os.path.exists(ledockinPath):
        shutil.rmtree(ledockinPath)
        os.mkdir(ledockinPath)
    else:
        os.mkdir(ledockinPath)

    receptorListFile = os.path.join(receptorpath, 'ledock_in.list')
    line_list = []
    if os.path.isfile(receptorListFile):
        ledock_file_open = open(receptorListFile, 'r')
        ledockfileopen = ledock_file_open.read()
        line_list = ledockfileopen.splitlines()
        ledock_file_open.close()
    else:
        receptorFileList = get_filelist(receptorpath)
        ledock_file_open = open(receptorListFile, 'a')
        for receptorFile in receptorFileList:
            if receptorFile.endswith('_ledock.in'):
                line_list.append(receptorFile)
                ledock_file_open.write(receptorFile + '\n')
        ledock_file_open.close()
    ledockConfigPathList = []
    try:
        for ledockConfigFile in line_list:
            # 确定中间的目录
            middlePath = os.path.dirname(ledockConfigFile.replace(receptorpath, ''))
            ledockinPathTemp = ledockinPath + middlePath
            os.makedirs(ledockinPathTemp)

            ledock_config_file_open = open(ledockConfigFile, 'r')
            ledock_config_file = ledock_config_file_open.read()
            line_list = ledock_config_file.splitlines()
            lines_num = len(line_list)
            ledock_config_file_open.close()

            for i in range(0, ligandNum):
                receptor_ledock_new = os.path.join(ledockinPathTemp, os.path.basename(ledockinPathTemp) + '_ledock' + str(i) + '.in')
                receptor_ledock_file = open(receptor_ledock_new, 'a')

                for line in range(0, lines_num):
                    if line == 1:
                        receptor_ledock_file.write(ledockConfigFile.replace('_ledock.in', '_pro.pdb') + '\n')
                    elif line == 15:
                        receptor_ledock_file.write(ligandFileListNew[i] + '\n')
                    else:
                        receptor_ledock_file.write(line_list[line] + '\n')
                receptor_ledock_file.close()
            ledockConfigPathList.append(ledockinPathTemp)

        used_time = time.time() - start_time
        logging.info("调整配置文件所使用时间：%s" % (str(datetime.timedelta(seconds=used_time))))

        return ledockConfigPathList
    except BaseException as e2:
        print(e2)


'''
分组配体文件
'''


def Split_file_ligand(num, currentPath, ligandFileList):
    start_time = time.time()

    # 配体文件列表目录
    ligandListPath = os.path.join(currentPath, 'ligandList')
    if os.path.exists(ligandListPath):
        shutil.rmtree(ligandListPath)
        os.mkdir(ligandListPath)
    else:
        os.mkdir(ligandListPath)

    lines_num = len(ligandFileList)

    ligand_count = int(lines_num / num)
    ligand_count_last = int(lines_num % num)

    ligandFileListNew = []

    for i in range(0, num):
        file_ligandpath = os.path.join(ligandListPath, 'ligandlist' + str(i))  # 创建新的ligandlist
        ligandFileListNew.append(file_ligandpath)
        file = open(file_ligandpath, 'a')
        for j in range(i * ligand_count, (i + 1) * ligand_count):
            getline = ligandFileList[j]
            file.write(getline + '\n')
        file.close()

    # TODO 余数放在最后一个文件
    if (ligand_count_last != 0):
        file = open(file_ligandpath, 'a')
        start_line = num * ligand_count
        end_line = lines_num
        for j in range(start_line, end_line):
            getline = ligandFileList[j]
            file.write(getline + '\n')
        file.close()

    used_time = time.time() - start_time

    logging.info("分组配体文件所使用时间：%s" % (str(datetime.timedelta(seconds=used_time))))
    return ligandFileListNew


'''
获取目录下所有文件
'''


def get_filelist(dir):
    Filelist = []
    for home, dirs, files in os.walk(dir):
        for filename in files:
            # 文件名列表，包含完整路径
            if filename == '.DS_Store':
                pass
            elif filename == 'ligandFile':
                pass
            else:
                Filelist.append(os.path.join(home, filename))
            # # 文件名列表，只包含文件名
            # Filelist.append( filename)

    return Filelist


'''
生成相关的文件
	截取配体文件路径中的最后一级目录名称作为临时目录
	根据配体文件列表，生成配体列表文件,多行进行拆分
	复制一份受体文件，配体文件到新生成的目录
	按照规则生成配置文件
	生成qsub脚本
		计算节点数
		计算核心数
		指定结果目录
		生成mpiexec条目
		所有运算结束之后，运行结果文件提取Python
'''


def main():
    logging.info('=' * 100)
    logging.info(str((datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))))
    logging.info('=' * 100)

    # 计算资源
    # 节点数
    nodes = 1
    # 核心数
    ppn = 36
    # 运算进程数
    running = nodes * ppn

    # 一个脚本文件中放多少个mpiexec任务
    # TODO 修改的时候，请确认一次能否提交这么多的任务
    # TODO 如果任务中包含的数量过多，服务器将拒绝提交该任务
    qsub_count = 1

    # TODO 测试用地址
    # currentPath = '/Volumes/data/home/Downloads/temp/task0001'

    # TODO 正式地址
    currentPath = os.getcwd()
    logging.info('工作目录：' + currentPath)

    # TODO 测试用地址
    # receptorpath = '/Volumes/data/home/Developer/workspace/centling/code/deep/ledock/Mpro/Mpro'

    # TODO 正式地址 高研院
    receptorpath = '/lustre/home/weizhiqiang/dev/zyd/20210420/Mpro'

    # TODO 正式地址 即墨超算
    # receptorpath = '/home/csywz/deep/Mpro'

    logging.info('受体文件目录：' + receptorpath)

    ligandPath = os.path.join(currentPath, 'ligand')

    # 对上传上来的配体压缩包解压
    ligandZipFile = os.path.join(ligandPath, 'ligandFile.zip')
    zipFile = zipfile.ZipFile(ligandZipFile)
    zipFile.extractall(ligandPath)
    zipFile.close()
    os.remove(ligandZipFile)

    # 生成 ligandList 文件列表
    ligandFileList = get_filelist(ligandPath)

    # 根据配体的数量来动态调整资源
    ligandSize = len(ligandFileList)

    # 上传上来的配体文件列表
    ligandFile = os.path.join(ligandPath, 'ligandFile')

    ligand_file_open = open(ligandFile, 'r')
    ligandfileopen = ligand_file_open.read()
    ligand_file_list = ligandfileopen.splitlines()
    ligand_file_open.close()
    ligand_file_open_length = len(ligand_file_list)

    if ((ligandSize) != ligand_file_open_length):
        message = '上传到服务器的配体文件数量和配体文件列表中的数量不相等'
        logging.error(message)
        return message

    # 校验一下数据有没有问题
    for ligand_file in ligand_file_list:
        ligand_file_name = os.path.join(ligandPath, os.path.basename(ligand_file))
        if (ligand_file_name not in ligandFileList):
            error_message = '本地文件：' + ligand_file + ' 在服务器上：' + ligandPath + '不存在'
            logging.error(error_message)
            return error_message

    #logging.info("配体文件校验通过：需要计算的配体文件数量：" + str(ligandSize))
    #if (ligandSize < running):
    #    running = 1
    #    nodes = 1
    #    ppn = 1

    # 生成切割后的配体文件列表
    ligandFileListNew = Split_file_ligand(running, currentPath, ligandFileList)

    # 生成配置文件
    ledockConfigPathList = Split_file_pdb(currentPath, receptorpath, ligandFileListNew)

    # 脚本目录
    qsubPath = os.path.join(currentPath, 'qsub_pbs')
    if os.path.exists(qsubPath):
        shutil.rmtree(qsubPath)
        os.mkdir(qsubPath)
    else:
        os.mkdir(qsubPath)

    # 生成脚本列表
    qsubPathList = receptor_qsub(nodes, ppn, qsub_count, currentPath, ledockConfigPathList)

    # 执行 提交任务
    qstatList = []
    for qsub in qsubPathList:
        workPath=os.path.dirname(qsub)
        qsunName= os.path.basename(qsub)
        os.chdir(workPath)

        # TODO 高研院
        qsub_command = '/usr/bin/qsub ' + qsunName

        # TODO 即墨超算
        # qsub_command = '/home/SystemSoftware/tsce/torque6/bin/qsub ' + qsunName

        result = os.system(qsub_command)
        logging.info('提交任务：'+qsub)
        qstatList.append(result)
    return qstatList


if __name__ == "__main__":
    main()
