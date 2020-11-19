#!/bin/bash

# Written by Mrinmoy S. Roddur in Fall 2020

#file locations

base="$HOME/scratch/ASTRAL-Pro_vs_FastMulRFS/INDELible_Pipeline/"
data="$base/sim_data"
code="$base/sim_code"

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

# param: the location of the directory
# loops through each replicate, creates two directory - 100 and 500, stores sequences and estimated trees in there, 
# and also makes a copy of the trees in trees_est files

function estimate_trees {
    for repl in ${repls[@]}; do
        #trees.txt file contains all the gene trees. Unfortunately, INDELible works with one file containing all the trees
        rm -f $1/$repl/trees.txt
        for gtree in ${total_gtrees[@]}; do
            if [ -f $1/$repl/g_trees$gtree\.trees ]; then
                cat $1/$repl/g_trees$gtree\.trees >> $1/$repl/trees.txt
            fi
        done
        for seqlen in ${sqlens[@]}; do
            rm -rf $1/$repl/$seqlen
            mkdir $1/$repl/$seqlen

            # takes as input the number of gene trees and sequence lengths, and generates some parameter files.
            python $setindeliblepy \
                -n $n_total_gtrees \
                -f 113.48869 69.02545 78.66144 99.83793 \
                -r 12.776722 20.869581 5.647810 9.863668 30.679899 3.199725 \
                -a -0.470703916 0.348667224 \
                -l $seqlen \
                -o $1/$repl/indelible-parameters.csv


            # based on the parameter file, generates all phylip files and store them in 100/ or 500/ directory
            python $indeliblepy \
                -x $indelible \
                -s 1 \
                -e $n_total_gtrees \
                -p $1/$repl/indelible-parameters.csv \
                -t $1/$repl/trees.txt \
                -o $1/$repl/$seqlen

	        rm -f $1/$repl/trees_est\_$seqlen.txt

            # runs FastTree on phylip files. I've also copied the trees in trees_est_$seqlen.txt files. 
            for gtree in ${total_gtrees[@]}; do
                $fasttree -nt -gtr $1/$repl/$seqlen/$gtree\.phy > $1/$repl/$seqlen/$gtree\.nwk
                $1/$repl/$seqlen/$gtree\.nwk >> $1/$repl/trees_est\_$seqlen.txt
            done
            rm $1/$repl/indelible-parameters.csv
        done
    done
}

estimate_trees $data/default

: '
for gdr in ${gdrs[@]}; do
    for glr in ${glrs[@]}; do
        if [ -d $data/exp_gdl_$gdr\_$glr ]; then
            estimate_trees $data/exp_gdl_$gdr\_$glr
        fi
    done
done

for ils in ${ilss[@]}; do
    estimate_trees $data/exp_ils_$ilss
done

estimate_trees $data/exp_species
'
