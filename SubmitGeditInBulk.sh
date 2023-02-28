#!/usr/bin/bash

# This script can be used to peform deconvolution of one bulk mix file with several reference matrixes, on the cluster
 # $1 path to mix file
 # $2 path to reference file or files
 # $3 path to directory to save the results


#MANUAL ARGUMENTS: these will be used if there are unused arguments when the script is called in the terminal
manual1='/omics/groups/OE0436/data/koppany/TCGA_processed/tpm_with_gene_names/*.tsv'
manual2='/omics/groups/OE0436/data/koppany/TISCH/reference_matrices/*.txt'
manual3='/home/g818y/koppany_link/deconv_results/TCGA_mix_TISCH_ref'

# START OF SCRIPT
#assign argument1
if [[ "$1" != "" ]]; then
    input1="$1"
else
    input1=$manual1
fi
#assign argument2
if [[ "$2" != "" ]]; then
    input2="$2"
else
    input2=$manual2
fi
#assign argument3
if [[ "$3" != "" ]]; then
    input3="$3"
else
    input3=$manual3
fi

#print argumnets
echo "TEST: arg1 = $input1"
echo "TEST: arg2 = $input2"
echo "TEST: arg3 = $input3"

#path to mix/sample files
mixArray=($(ls $input1))
echo "mix/sample files to be deconvoluted: "
for name in "${mixArray[@]}" 
do
echo "$name"
done

#path to reference matrixes.
refArray=($(ls $input2))
echo "references to be used:"

for name in "${refArray[@]}" 
do
echo "$name"
done

# the directory where the deconvolution results should be saved
outDir=$input3 
echo "Directory where the files will be saved: $outDir"

if [ -d "$outDir" ]; then
  ### Take action if $outDIR exists ###
  echo "placing deconvolution results in already existing $outDir"
else
  ###  Control will jump here if $DIR does NOT exists ###
  echo " ${outDir} will be created"
  mkdir $outDir
fi

   #intoduce counter to submit only up to 200 jobs at a time
   nsub=0

for currentMix in "${mixArray[@]}"
    do
    # get the name of the mix by removing the extension and dir path
    mixname=${currentMix##*/}
    mixname=${mixname%.*}

    for currentRef in "${refArray[@]}"
        do 
        #get the name of the reference, by removing the extension and directry name
        refname=${currentRef##*/}
        refname=${refname%.*}
        outname=${outDir}/${mixname}_${refname}.tsv
        jobname=${mixname}_${refname}_GEDIT_deconv
    
        #perform deconvolution by GEDIT
        echo "current Mix: $mixname, current reference: $refname"
        echo "Result will be saved as: $outname" 
        echo "Submitting $jobname to cluster"
        
         nsub=$(($nsub+1)) # keep tranck of how many jobs were submitted
        echo "nsub: $nsub "
        #submit job
        bsub -J $jobname \
        -W 0:20 \
        -o "/home/g818y/koppany_link/lsf_log/${jobname}_out.log" \
        -e "/home/g818y/koppany_link/lsf_log/${jobname}_error.log" \
        -n 1 -R "rusage[mem=16G]" \
        python GEDIT2.py -mix $currentMix -ref $currentRef -outFile $outname 
        
        #check submissions in every 200th occasion
        if (( $nsub > 200 )); then
            echo '200th job was submitted. waiting for 5 minutes '
            sleep 5m 
            echo 'Currently active jobs:'
            bjobs
            
            echo 'Number of active jobs:'
            bjobs | wc -l
            
            echo 'continue? [y/n]:'
            read ans1
            
            while [ $ans1 == n ] ; do
                echo 'sleeping for 3 minutes'
                sleep 3m
                
                echo 'Currently active jobs:'
                bjobs
                echo 'Number of active jobs:'
                bjobs | wc -l
                echo 'continue? [y/n]:'
                read ans1
            done 
            echo 'submitting the next 200 jobs'
            nsub=0 #restart nsub counter
        fi
        
    done
done
