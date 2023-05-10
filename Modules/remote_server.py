#!/usr/bin/python
# coding=utf-8
import datetime
import os
import time

import logging
import zipfile

logging.basicConfig(
    filename='remote_server.log',
    level=logging.DEBUG,
    format=' %(asctime)s - %(levelname)s - [文件：%(filename)s 第%(lineno)d行 方法：%(funcName)s] - %(message)s'
)
logging.captureWarnings(True)

import paramiko


# 校验目前队列状态
def checkJobStat(sshClient, jobId):
    # TODO 高研院
    checkJobStatCommand = "/usr/bin/qstat " + jobId

    # TODO 即墨超算
    # checkJobStatCommand = "/home/SystemSoftware/tsce/torque6/bin/qstat " + jobId

    ststus = True
    try:
        stdin, stdout, stderr = sshClient.exec_command(checkJobStatCommand);
        result_stdout = stdout.read()
        result_stdout_str = result_stdout.decode('utf-8')
        if len(result_stdout_str) > 0:
            print(result_stdout.decode('utf-8'))
            if ' R ' in str(result_stdout):
                print('检测到 %s 的运行状态 为：  R ' % (jobId))
            elif ' E ' in str(result_stdout):
                print('检测到 %s 的运行状态 为：  E  ' % (jobId))
            elif ' C ' in str(result_stdout):
                ststus = False
                print('检测到 %s 的运行状态 为：  C ' % (jobId))
            elif ' Q ' in str(result_stdout):
                print('检测到 %s 的运行状态 为：  Q  ' % (jobId))
            else:
                print('检测到 %s 的运行状态为其他！' % (jobId))
                ststus = False
        else:
            print('没有检测到 %s 的运行状态' % (jobId))
            ststus = False
    except BaseException as e1:
        print(e1)
    finally:
        return ststus


'''
==========================
=请调用该方法,并传入相关的参数=
==========================
'''


# _ligandListFilePath 本地配体文件列表
# _remotepath 服务器有可写权限的目录

def sftpClient(
        _ligandListFilePath,
        _remotepath='/home/test/',
        _hostname='10.211.55.5',
        _port=22,
        _username='test',
        _passwd='1234.abcd'
):
    _start_time = time.time()
    logging.info('=' * 100)
    logging.info(str((datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))))
    logging.info('=' * 100)

    logging.info('本地配体列表文件:' + _ligandListFilePath)


    # 工作空间(taskooo1)
    remoteServerTaskPathName = os.path.basename(os.path.dirname(_ligandListFilePath))

    # 服务器工作目录
    remoteServerTaskPath = _remotepath + remoteServerTaskPathName
    logging.info('服务器工作目录:' + remoteServerTaskPath)

    # 服务器配体文件路径
    remoteLigandPath = os.path.join(remoteServerTaskPath, 'ligand')
    logging.info('服务器配体文件路径:' + remoteLigandPath)

    currentPath = os.getcwd()+'/Modules'
    #currentPath = os.getcwd()

    # 任务列表
    localTaskPath = _ligandListFilePath.rpartition("/")[0]
    subsFile = os.path.join(localTaskPath, 'subs')
    getResultsubsFile = os.path.join(localTaskPath, 'getResultsubs')
    remoteResultFile = os.path.join(remoteServerTaskPath, 'result', 'result.md')
    localResultFile = os.path.join(localTaskPath, 'result.md')
    try:
        # SSH 客户端
        sshClient = paramiko.SSHClient()
        sshClient.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        sshClient.connect(hostname=_hostname, port=_port, username=_username, password=_passwd, compress=True)

        # SFTP 客户端
        sftp = sshClient.open_sftp()

        print('正在连接主机%s......' % (_hostname))


    except BaseException as e:
        print(e)
        print("服务器链接错误！请检查相关配置")
    else:
        print('成功连接主机%s......' % (_hostname))

        try:
            if not os.path.exists(getResultsubsFile):
                if not os.path.exists(subsFile):
                    # 压缩配体文件
                    zipligandListFilePath = _ligandListFilePath + '.zip'
                    zipFile = zipfile.ZipFile(zipligandListFilePath, 'w', zipfile.ZIP_DEFLATED)
                    ligandListNewPath = os.path.join(os.path.dirname(_ligandListFilePath), 'ligandFile')
                    # 本地配体文件重命名
                    os.rename(_ligandListFilePath, ligandListNewPath)
                    zipFile.write(ligandListNewPath, os.path.basename(ligandListNewPath))

                    ligandListPathFile = open(ligandListNewPath, 'r')
                    ligandListPathFileLines = ligandListPathFile.readlines()
                    ligandListPathFile.close()

                    for line in ligandListPathFileLines:
                        rs = line.strip()
                        print('正在压缩文件：', rs)
                        zipFile.write(rs, os.path.basename(rs))
                    zipFile.close()
                    time.sleep(4)

                    os.rename(ligandListNewPath, _ligandListFilePath)
                    sshClient.exec_command("mkdir -vp " + remoteLigandPath)
                    time.sleep(3)  # 创建文件夹需要时间
                    start_upload_time = time.time()

                    print('正在向主机' + _hostname + '上上传配体列表文件: ' + zipligandListFilePath + ' ......')
                    sftp.put(zipligandListFilePath,os.path.join(remoteLigandPath, 'ligandFile.zip'))
                    logging.info('上传：' + zipligandListFilePath)

                    end_upload_time = time.time()

                    uploadUsedTime = end_upload_time - start_upload_time
                    uploadTimeStr = '上传配体文件使用的总时间：' + str(datetime.timedelta(seconds=uploadUsedTime))
                    print(uploadTimeStr)
                    logging.info(uploadTimeStr)

                    # 上传计算脚本
                    sftp.put(os.path.join(currentPath, 'prepare_files.py'), os.path.join(remoteServerTaskPath, 'prepare_files.py'))
                    print('正在向主机' + _hostname + '上上传 Python 文件: ' + os.path.join(currentPath, 'prepare_files.py ......'))

                    sftp.put(os.path.join(currentPath, 'mpi_ledock.py'), os.path.join(remoteServerTaskPath, 'mpi_ledock.py'))
                    print('正在向主机' + _hostname + '上上传 Python 文件: ' + os.path.join(currentPath, 'mpi_ledock.py ......'))

                    sftp.put(os.path.join(os.path.dirname(currentPath), 'ledock', 'dock.sh'), os.path.join(remoteServerTaskPath, 'dock.sh'))
                    print('正在向主机' + _hostname + '上上传 dock.sh 文件: ' + os.path.join(currentPath, 'dock.sh ......'))
                    sshClient.exec_command('chmod +x ' + os.path.join(remoteServerTaskPath, 'dock.sh'))

                    sftp.put(os.path.join(currentPath, 'delete.sh'), os.path.join(remoteServerTaskPath, 'delete.sh'))
                    print('正在向主机' + _hostname + '上上传 delete.sh 文件: ' + os.path.join(currentPath, 'delete.sh ......'))
                    sshClient.exec_command('chmod +x ' + os.path.join(remoteServerTaskPath, 'delete.sh'))

                    # 生成执行脚本
                    localRunBush = os.path.join(currentPath, 'run.sh')
                    remoteRunBush = os.path.join(remoteServerTaskPath, 'run.sh')

                    file_runBush = open(localRunBush, 'w')
                    file_runBush.write('#!/bin/bash' + '\n')

                    file_runBush.write('cd ' + remoteServerTaskPath + '\n')
                    file_runBush.write('python prepare_files.py' + '\n')
                    file_runBush.close()

                    sftp.put(localRunBush, remoteRunBush)
                    print('正在向主机' + _hostname + '上上传 run.sh 文件: ' + localRunBush + '......')
                    time.sleep(4)
                    sshClient.exec_command('chmod +x ' + remoteRunBush)

                    print('在主机' + _hostname + '上执行: ' + remoteRunBush + '......')

                    isNotWork = True
                    while isNotWork:
                        time.sleep(3)  # 等待三秒之后执行
                        stdin, stdout, sterr = sshClient.exec_command(remoteRunBush)
                        # 检测运行状态
                        qsubList = stdout.readlines()
                        if qsubList:
                            if '.mu01' in qsubList[1]:
                                isNotWork = False
                        else:
                            isNotWork = True

                    subsFileOpen = open(subsFile, 'w')
                    subsFileOpen.writelines(qsubList)
                    subsFileOpen.close()
                else:
                    subsFileOpen = open(subsFile, 'r')
                    qsubList = subsFileOpen.readlines()
                    subsFileOpen.close()
                    logging.info("存在日志文件！")

                start_time = time.time()
                for qsubName in qsubList:
                    if ":" in qsubName:
                        continue
                    list_start = qsubList.index(qsubName)
                    qsubList = qsubList[list_start:]
                    qsubName = qsubName.strip()
                    while checkJobStat(sshClient, qsubName):
                        time.sleep(2)

                end_time = time.time()
                print('对接运行时间为： %s' % (str(datetime.timedelta(seconds=(end_time - start_time)))))
                # 获取结果
                # 上传获取结果脚本
                sftp.put(os.path.join(currentPath, 'get_result.py'), os.path.join(remoteServerTaskPath, 'get_result.py'))
                print('正在向主机' + _hostname + '上上传 Python 文件: ' + os.path.join(currentPath, 'get_result.py ......'))

                sftp.put(os.path.join(currentPath, 'mpi_getresult.py'), os.path.join(remoteServerTaskPath, 'mpi_getresult.py'))
                print('正在向主机' + _hostname + '上上传 Python 文件: ' + os.path.join(currentPath, 'mpi_getresult.py ......'))

                # 生成获取结果执行脚本
                localRunBush = os.path.join(currentPath, 'runGetResult.sh')
                remoteRunBush = os.path.join(remoteServerTaskPath, 'runGetResult.sh')

                file_runBush = open(localRunBush, 'w')
                file_runBush.write('#!/bin/bash' + '\n')
                file_runBush.write('cd ' + remoteServerTaskPath + '\n')
                file_runBush.write('python get_result.py ' + os.path.join(remoteServerTaskPath, 'result') + '\n')
                file_runBush.close()

                sftp.put(localRunBush, remoteRunBush)
                print('正在向主机' + _hostname + '上上传 runGetResult.sh 文件: ' + localRunBush + '......')
                time.sleep(4)
                sshClient.exec_command('chmod +x ' + remoteRunBush)

                print('在主机' + _hostname + '上执行: ' + remoteRunBush + '......')

                # TODO 获取生成之后的结果
                isNotWork = True
                while isNotWork:
                    time.sleep(3)  # 等待三秒之后执行
                    stdin, stdout, sterr = sshClient.exec_command(remoteRunBush)
                    qsubList = stdout.readlines()
                    if qsubList:
                        if '.mu01' in qsubList[0]:
                            isNotWork = False
                    else:
                        isNotWork = True

                # 检测运行状态
                getResultsubsFileOpen = open(getResultsubsFile, 'w')
                getResultsubsFileOpen.writelines(qsubList)
                getResultsubsFileOpen.close()

            else:
                getResultsubsFileOpen = open(getResultsubsFile, 'r')
                qsubList = getResultsubsFileOpen.readlines()
                getResultsubsFileOpen.close()
            _get_result_start_time = time.time()
            # 检测获取结果的任务
            for qsubName in qsubList:
                qsubName = qsubName.strip()
                while checkJobStat(sshClient, qsubName):
                    time.sleep(2)

            print('正在下载结果文件到本地......')
            time.sleep(4)
            sftp.get(remoteResultFile, localResultFile)

            resultFile = open(localResultFile, 'r')
            result_list_str = resultFile.read()
            resultFile.close()

            _get_result_end_time = time.time()

            _get_result_time = _get_result_end_time - _get_result_start_time
            print('获取结果使用的时间：%s' % (str(datetime.timedelta(seconds=_get_result_time))))

            _end_time = time.time()

            _test_time = _end_time - _start_time

            print('计算用总时间：%s' % (str(datetime.timedelta(seconds=_test_time))))
            print(result_list_str)
            print(len(result_list_str))
            return result_list_str

        except BaseException as e1:
            print(e1)
    finally:
        sshClient.close()
        sftp.close()


# 配体文件列表，其中 task0001 为任务编号
# /Volumes/data/home/Developer/workspace/centling/code/deep/ledock/ligandlist/task0001/ligand02

def main(ligandListFilePath):
    print('=' * 67)
    print('=' * 10, '该方法仅供测试使用！请调整好相关参数,请不要用于生产环境', '=' * 10)
    print('=' * 67)
    # 高研院
    sftpClient(
        _ligandListFilePath=ligandListFilePath,
        _remotepath='/lustre/home/weizhiqiang/',
        _hostname='10.130.2.222',
        _username='weizhiqiang',
        _passwd='IAOS1234'
    )

    # 即墨超算
    # sftpClient(
    #     _ligandListFilePath=ligandListFilePath,
    #     _remotepath='/home/csywz/',
    #     _hostname='11.11.100.10',
    #     _username='csywz',
    #     _passwd='csywz@20131218'
    # )

    # sftpClient(
    #     _ligandListFilePath=ligandListFilePath
    # )


# 调用测试
if __name__ == '__main__':
    print('=' * 35)
    print('=' * 10, '开始运行测试程序', '=' * 10)
    print('=' * 35)
    ligandListFilePath = '/Users/macbookpro/rldmpo/ledock/ligandlist/task0001/ligandFile'
    main(ligandListFilePath)
