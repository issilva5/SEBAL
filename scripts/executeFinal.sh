for i in $(seq 3 30);
do
    sh execute.sh &&	

    for j in $(ls /home/ubuntu/TDir | grep r_tmp);
    do
	echo "Removing " ${j}
	rm /home/ubuntu/TDir/${j}
    done

    echo "Compactando..."
    tar -czvf teste${i}.tar.gz /home/ubuntu/TDir/ &&

    echo "Movendo para diretorio final..."
    mv teste${i}.tar.gz /home/ubuntu/Output/ &&

    for j in $(ls /home/ubuntu/TDir | grep .nc);
    do
	echo "Removing " ${j}
	rm /home/ubuntu/TDir/${j}
    done

    echo "Compactando apenas arquivos csv..."
    tar -czvf testedocs${i}.tar.gz /home/ubuntu/TDir/ &&

    echo "Movendo para diretorio final..."
    mv testedocs${i}.tar.gz /home/ubuntu/Output/ &&

    for j in $(ls /home/ubuntu/TDir);
    do
	echo "Removing " ${j}
	rm /home/ubuntu/TDir/${j}
    done

    echo "Finalizando teste " ${i}
done
