#!/bin/bash
for fn in /Volumes/EBSExtDrive/RNAseq/UNC130_Sep7-2_S3_L005;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i elegans_index -l A \
         -1 ${fn}/${samp}_01.fastq \
         -2 ${fn}/${samp}_02.fastq \
		 --gcBias True -o /Users/emilybeckettsward/Desktop/emilysalmonrun/quants/${samp}_quant
done 
