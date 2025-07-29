#!/bin/bash

indir=$1
outfile=$2
cd $indir
find ./ -name "*report*" -type f -exec tar rf $outfile {} \;
tar rf $outfile msa/*