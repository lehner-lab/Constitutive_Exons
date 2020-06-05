#!/bin/bash
#Export all enviroment variables
#$ -V
#Use current working directory
#$ -cwd
#Join stdout and stderr
#$ -j y
#Email address
#$ -m ae
#$ -M pablo.baeza@crg.eu
#Memory
#$ -l virtual_free=10G
#Time
#$ -q short-sl7,long-sl7
#$ -l h_rt=6:00:00
#Parallel environment
#$ -pe smp 4
#$ -N dimsum_FAS_Human
#$ -t 1-1

fastqFileDir="full path to directory with the dataset"
fastqFileExtension=".fastq"
gzipped="TRUE"
experimentDesignPath="full path to file named experimentDesign_FAS_Human.txt"
barcodeDesignPath="full path to file named barcodeDesign_FAS_Human.txt"
barcodeErrorRate="0"
cutadaptMinLength="50"
cutadaptErrorRate="0.2"
usearchMinQual="30"
usearchMaxee="0.5"
outputPath="full path to directory where DiMSum will output the processed files"
projectName="FAS_Human"
startStage="0" # if the data is already demultiplexed into replicates, set this to "1"
wildtypeSequence="GATCCAGATCTAACTTGGGGTGGCTTTGTCTTCTTCTTTTGCCAATTCCACTAATTGTTTGGG"
numCores="4"
usearchMinlen="50"
stranded="F"
experimentDesignPairDuplicates="T"
retainIntermediateFiles="T"
cutadaptCut3First="5"
cutadaptCut3Second="5"
fitnessMinInputCountAll="3000"

DiMSum -i $fastqFileDir -l $fastqFileExtension -g $gzipped -e $experimentDesignPath -b $barcodeDesignPath -n $cutadaptMinLength -q $usearchMinQual -m $usearchMaxee -o $outputPath -p $projectName -s $startStage -w $wildtypeSequence -c $numCores --usearchMinlen $usearchMinlen --stranded $stranded --barcodeErrorRate $barcodeErrorRate --experimentDesignPairDuplicates $experimentDesignPairDuplicates --retainIntermediateFiles $retainIntermediateFiles --cutadaptCut3First $cutadaptCut3First --cutadaptCut3Second $cutadaptCut3Second --fitnessMinInputCountAll $fitnessMinInputCountAll

