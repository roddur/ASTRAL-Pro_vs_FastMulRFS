#!/bin/bash

# Written by Mrinmoy S. Roddur in Fall 2020

#File locations

base="$HOME/ASTRAL-Pro_vs_FastMulRFS/INDELible_Pipeline"
data="$base/sim_data"
code="$base/sim_code"
batchfile="$code/template.sbatch"

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

function estimate_trees_parallel {

    # Takes as input the name of the dataset directory

    for repl in ${repls[@]}; do
        rm -f $data/$1/$repl/ind_job.sbatch
        touch $data/$1/$repl/ind_job.sbatch
        cat $batchfile > $data/$1/$repl/ind_job.sbatch
        sed -i "s/TEMPLATE/$1\_$repl\_job/g" $data/$1/$repl/ind_job.sbatch
        echo $code/job_per_repl.sh $1 $repl >> $data/$1/$repl/ind_job.sbatch
        sbatch $data/$1/$repl/ind_job.sbatch
    done

}

estimate_trees_parallel $1
