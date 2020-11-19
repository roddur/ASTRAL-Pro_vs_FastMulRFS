#!/bin/bash

# Written by Mrinmoy S. Roddur in Fall 2020

#File locations

base="$HOME/Projects/Compare_FastMulRFS_with_ASTRAL-Pro/INDELible_Pipeline"
data="$base/sim_data"
code="$base/sim_code"
batchfile="job.sbatch"

setindeliblepy="$code/set_indelible_params.py"
indeliblepy="$code/run_indelible.py"
indelible="$code/indelible"
fasttree="$code/FastTree"

#Constants

nrepl=10
repls=( $( seq -f "%02g" 1 $nrepl ) )
n_total_gtrees=1000
total_gtrees=( $( seq -f "%04g" 1 $n_total_gtrees ) )
gdrs=( 0.000000001 0.0000000001  0.0000000005 )
glrs=( 0 05 1 )
ilss=( 0 2e8 5e8)
sqlens=( 100 500 )
n_gtrees=( 100 1000)

function job_per_replicate {
    
    # generates all estimated trees for a specific replicate
    # param the location of the replicate

    rm -f $data/$1/trees.txt
    for gtree in ${total_gtrees[@]}; do
        if [ -f $data/$1/g_trees$gtree\.trees ]; then
            cat $data/$1/g_trees$gtree\.trees >> $1/$repl/trees.txt
        fi
    done

    for seqlen in ${sqlens[@]}; do
        rm -rf $data/$1/$seqlen
        mkdir $data/$1/$seqlen

        # takes as input the number of gene trees and sequence lengths, and generates some parameter files.
        echo python $setindeliblepy \
            -n $n_total_gtrees \
            -f 113.48869 69.02545 78.66144 99.83793 \
            -r 12.776722 20.869581 5.647810 9.863668 30.679899 3.199725 \
            -a -0.470703916 0.348667224 \
            -l $seqlen \
            -o $data/$1/indelible-parameters.csv


        # based on the parameter file, generates all phylip files and store them in 100/ or 500/ directory
        echo python $indeliblepy \
            -x $indelible \
            -s 1 \
            -e $n_total_gtrees \
            -p $data/$1/indelible-parameters.csv \
            -t $data/$1/trees.txt \
            -o $data/$1/$seqlen

	    rm -f $data/$1/$repl/trees_est\_$seqlen.txt

        # runs FastTree on phylip files. I've also copied the trees in trees_est_$seqlen.txt files. 
        for gtree in ${total_gtrees[@]}; do
            $fasttree -nt -gtr $data/$1/$seqlen/$gtree\.phy > $data/$1/$seqlen/$gtree\.nwk
            $data/$1/$seqlen/$gtree\.nwk >> $data/$1/trees_est\_$seqlen.txt
        done
        rm -f $data/$1/indelible-parameters.csv
    done

}

job_per_replicate $1