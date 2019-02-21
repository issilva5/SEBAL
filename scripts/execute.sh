sh executeRScript.sh ../workspace/R/Algoritmo-completo-v24012019.R /home/ubuntu/Dir/ /home/ubuntu/TDir/ &
sh collect-cpu-usage.sh $(pidof R) > /home/ubuntu/TDir/cpu.csv &
sh collect-memory-usage.sh $(pidof R) > /home/ubuntu/TDir/mem.csv &