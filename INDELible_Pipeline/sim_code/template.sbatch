#!/bin/bash
#
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --job-name=TEMPLATE
#SBATCH --partition=secondary
#SBATCH --output=out/TEMPLATE.o%j
##SBATCH --error=error/TEMPLATE.e%j
##SBATCH --mail-user=mroddur2@illinois.edu
##SBATCH --mail-type=BEGIN,END
#
# End of embedded SBATCH options
#