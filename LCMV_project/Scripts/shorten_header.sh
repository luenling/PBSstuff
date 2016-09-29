#!/bin/bash
#----------
# author: Lukas Endler
# date: 17.9.2016 at 16:34 
# Time-stamp: <2016-03-15 10:33:43 lukasendler>
# gets a
# call with bash command blub_1.fq blub_2.fq >> logfile.log 2>> log.error.log
# also sort the file afterwards
#--------------
BASEDIR=/Volumes/Temp/Lukas/LCMV_project
REFGENOME=$BASEDIR/References/viruses_Mus_musculus.GRCm38.fa.gz
PICARD=/usr/local/Cellar/picard-tools/1.128/share/java/picard.jar
SAMTOOLS=/usr/local/bin/samtools


