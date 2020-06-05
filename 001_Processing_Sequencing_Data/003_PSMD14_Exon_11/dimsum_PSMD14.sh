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
#$ -l virtual_free=20G
#Time
#$ -q long-sl7
#$ -l h_rt=99:00:00
#Parallel environment
#$ -pe smp 4
#$ -N dimsum_PSMD14
#$ -t 1-1

fastqFileDir="full path to directory with the dataset"
fastqFileExtension=".fastq"
gzipped="TRUE"
experimentDesignPath="full path to file named experimentDesign_PSMD14.txt"
barcodeDesignPath="full path to file named barcodeDesign_PSMD14.txt"
barcodeErrorRate="0"
cutadaptMinLength="50"
cutadaptErrorRate="0.2"
usearchMinQual="30"
usearchMaxee="0.5"
outputPath="full path to directory where DiMSum will output the processed files"
projectName="PSMD14_exon_11"
startStage="0"
wildtypeSequence="GCTGTAGAAGAAGAAGATAAGATGACACCTGAACAGCTGGCAATAAAGAATGTTGGCAAGCAG"
numCores="4"
usearchMinlen="63"
stranded="F"
experimentDesignPairDuplicates="T"
retainIntermediateFiles="T"
sequenceType="noncoding"


DiMSum -i $fastqFileDir -l $fastqFileExtension -g $gzipped -e $experimentDesignPath -b $barcodeDesignPath -n $cutadaptMinLength -q $usearchMinQual -m $usearchMaxee -o $outputPath -p $projectName -s $startStage -w $wildtypeSequence -c $numCores --usearchMinlen $usearchMinlen --stranded $stranded --experimentDesignPairDuplicates $experimentDesignPairDuplicates --barcodeErrorRate $barcodeErrorRate --retainIntermediateFiles $retainIntermediateFiles --sequenceType $sequenceType

