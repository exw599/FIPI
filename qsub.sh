#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 4
#$ -l h_vmem=0.5G
#$ -l h_rt=48:00:00 
#$ -m be
#$ -M cgu599@outlook.com
./run/AAA
