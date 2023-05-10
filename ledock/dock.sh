#!/bin/sh
# $1 结果文件夹
# $2 ledock 配置文件
dir=$1
if [ ! -d  $dir ]
then
	mkdir -p $dir
fi
cd $dir

# 高研院
/lustre/home/weizhiqiang/ledock/QMarine/ledock_go $2

# 即墨
# /home/csywz/deep/ledock/ledock_go $2
