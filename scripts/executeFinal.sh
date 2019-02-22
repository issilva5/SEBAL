#!/bin/bash

for i in `seq 2 2`;
do
    sh execute.sh
    rm /home/ubuntu/TDir/r*
    tar -czvf teste$i.tar.gz /home/ubuntu/TDir/
    mv teste$i.tar.gz /home/ubuntu/output/
    rm /home/ubuntu/TDir/*.*
done