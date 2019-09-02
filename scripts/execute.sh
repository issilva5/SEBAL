#!/bin/bash

sh executeRScript.sh ../workspace/R/Algoritmo-completo-v24012019.R /home/ubuntu/Dir/ /home/ubuntu/TDir/ &

sleep 1

PID=$(pidof R)
echo ${PID} > /home/ubuntu/pid

sh collect-cpu-usage.sh ${PID} > /home/ubuntu/TDir/cpu.csv &
sh collect-memory-usage.sh ${PID} > /home/ubuntu/TDir/mem.csv &
sh collect-disk-usage.sh ${PID} > /home/ubuntu/TDir/disk.csv
