#!/bin/bash

# Written by Mrinmoy S. Roddur in Fall 2020

base="/home/mrinmoy/Projects/Compare_FastMulRFS_with_ASTRAL-Pro/INDELible_Pipeline"
data="$base/sim_data"
code="$base/sim_code"

setindeliblepy="$code/set_indelible_params.py"
indeliblepy="$code/run_indelible.py"
indelible="$code/indelible"
fasttree="$code/FastTree"

nrepl=10
repls=( $( seq -f "%02g" 1 $nrepl ) )
n_total_gtrees=1000
total_gtrees=( $( seq -f "%04g" 1 $n_total_gtrees ) )
gdrs=( 0.000000001 0.0000000001  0.0000000005 )
glrs=( 0 05 1 )
ilss=( 0 2e8 5e8)
sqlens=( 100 500 )
n_gtrees=( 100 1000)

function set_indelible_params {
    for repl in ${repls[@]}; do
        rm -f $1/$repl/trees.txt
        for gtree in ${total_gtrees[@]}; do
            if [ -f $1/$repl/g_trees$gtree\.trees ]; then
                cat $1/$repl/g_trees$gtree\.trees >> $1/$repl/trees.txt
            fi
        done
        for seqlen in ${sqlens[@]}; do
            rm -rf $1/$repl/$seqlen
            mkdir $1/$repl/$seqlen
            python $setindeliblepy \
                -n $n_total_gtrees \
                -f 113.48869 69.02545 78.66144 99.83793 \
                -r 12.776722 20.869581 5.647810 9.863668 30.679899 3.199725 \
                -a -0.470703916 0.348667224 \
                -l $seqlen \
                -o $1/$repl/indelible-parameters.csv

            python $indeliblepy \
                -x $indelible \
                -s 1 \
                -e $n_total_gtrees \
                -p $1/$repl/indelible-parameters.csv \
                -t $1/$repl/trees.txt \
                -o $1/$repl/$seqlen

            for gtree in ${total_gtrees[@]}; do
                $fasttree -nt -gtr $1/$repl/$seqlen/$gtree\.phy > $1/$repl/$seqlen/$gtree\.nwk
                $1/$repl/$seqlen/$gtree\.nwk >> $1/$repl/trees_est\_$seqlen.txt
            done
            rm $1/$repl/indelible-parameters.csv
        done
    done
}

set_indelible_params $data/default 100

for gdr in ${gdrs[@]}; do
    for glr in ${glrs[@]}; do
        if [ -d $data/exp_gdl_$gdr\_$glr ]; then
            set_indelible_params $data/exp_gdl_$gdr\_$glr 100
        fi
    done
done

for ils in ${ilss[@]}; do
    set_indelible_params $data/exp_ils_$ilss 100
done

set_indelible_params $data/exp_species 1000
