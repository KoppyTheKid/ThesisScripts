#!/usr/bin/bash

for ((i = 1 ; i < 3 ; i++)); do

jobname="test_${i}_job"

 bsub -J $jobname \
        -W 0:05 \
        -o "/home/g818y/koppany_link/lsf_log/${jobname}_out.log" \
        -e "/home/g818y/koppany_link/lsf_log/${jobname}_error.log" \
        -n 1 -R "rusage[mem=1G]" \
        sleep 10s
        
  if (( $i == 1 )); then
  
  echo 'Active jobs:'
  bjobs
    echo 'Number of active jobs:'
    bjobs | wc -l
    
  echo 'continue? [y/n]:'
  read ans1
  
   while [ $ans1 == n ] ; do
   echo 'sleeping for 10s'
   
   sleep 5s
   
    bjobs
    echo 'Number of active jobs:'
    bjobs | wc -l
    echo 'continue? [y/n]:'
    read ans1
   
   done   
  fi 
done
