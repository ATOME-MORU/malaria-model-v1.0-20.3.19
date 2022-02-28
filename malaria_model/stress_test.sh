#!/bin/bash

num_processes=

usage() {
	echo "usage: stress_test.sh [ [-p num_processes] | -h]"
}

while [ "$1" != "" ]; do
    case $1 in
        -p | --processes )      shift
                                num_processes=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done


for (( i = 0; i < $num_processes; i++ )); do
	device=$[ $i % 2 ]
	echo "starting $i on $device" 
	# sleep 1
	./build/bin/malaria_model_release_cuda \
		-c config.test.json \
		-o temp \
		-i run$i \
		-d $device &> ./temp/'run'$i'.log' &

    pids[${i}]=$!
done

for pid in ${pids[*]}; do
    wait $pid
done

# Once all processes are finished, run:
# grep "Time" temp/*.log