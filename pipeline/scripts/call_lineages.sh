#!/bin/bash
fa=$1
outfile=$2
logfile=$3
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pangolin
pangolin --update
pangolin $1 --outfile $2 > $3