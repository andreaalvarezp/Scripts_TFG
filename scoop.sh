#!/bin/bash
#
#SBATCH --job-name=ssp
#SBATCH --mem=100GB
#SBATCH -n 20
#SBATCH -p test
#SBATCH --ntasks=1
#SBATCH --time=10:00

## Defino variables
db="/home/aalvarez/SSP_final/db"
query="/home/aalvarez/SCOOP_final"
results="/home/aalvarez/SCOOP_final/results"
temp="/home/aalvarez/SCOOP_final/temp"
tables="/home/aalvarez/SCOOP_final/tables"
db_ssp="/home/aalvarez/SSP_final/db"
query_ssp="/home/aalvarez/SSP_final/query"
ssp="/home/aalvarez/SSP_final"
results_ssp="/home/aalvarez/SSP_final/results"
temp_ssp="/home/aalvarez/SSP_final/temp"
dir=$PWD

## Cargo módulos

module load BLAST+/2.9.0-gompi-2019a

### PEPTIDOS SCOOP

echo "Búsqueda de motivos SCOOP en $1"
echo "Posibles argumentos: scoop10, scoop12, fusarium, verticilium, magnaporthe, neurospora, proscoop"

# 1. Creo base de datos con el archivo fasta del proteoma de las tres cepas de Plectosphaerella
makeblastdb -in $db/final_scaffolds.functionalproteins.fasta -out $db/NR/final_scaffolds_db -dbtype 'prot'

# 2. BLAST-P 
blastp -db $db/NR/final_scaffolds_db -query $query/$1.fasta -out $results/$1_align

# Extraigo archivo fasta

grep -e ">" $results/$1_align > $temp/$1-Scaffoldslistnew.txt
cat $temp/$1-Scaffoldslistnew.txt | sort | uniq > $temp/$1-Scaffoldslistnewuni.txt
cut  -d " " -f 1 $temp/$1-Scaffoldslistnewuni.txt > $temp/$1-list1.txt
awk '{if($0 ~ /^>/){print $0} else {printf $0}}' $db/final_scaffolds.functionalproteins.fasta | perl -pe "s/>/\n>/g" > $temp/$1_oneline.fasta
grep -A1 --no-group-separator -f $temp/$1-list1.txt $temp/$1_oneline.fasta > $tables/$1-significant.fasta
echo "Number in fasta file"
grep -o ">" $tables/$1-significant.fasta | wc -l

### PEPTIDOS SSP

echo "Búsqueda de motivos SSP en $1"
echo "Posibles argumentos: ssp1-ssp14"

# 1. Creo base de datos con el archivo fasta del proteoma de las tres cepas de Plectosphaerella
makeblastdb -in $db/final_scaffolds.functionalproteins.fasta -out $db/NR/final_scaffolds_db -dbtype 'prot'

# 2. BLAST-P 
blastp -db $db_ssp/NR/final_scaffolds_db -query $query_ssp/$1.fasta -out $results_ssp/$1_align

# Extraigo archivo fasta

grep -e ">" $results_ssp/$1_align > $temp_ssp/$1-Scaffoldslistnew.txt
cat $temp_ssp/$1-Scaffoldslistnew.txt | sort | uniq > $temp_ssp/$1-Scaffoldslistnewuni.txt
cut  -d " " -f 1 $temp_ssp/$1-Scaffoldslistnewuni.txt > $temp_ssp/$1-list1.txt
awk '{if($0 ~ /^>/){print $0} else {printf $0}}' $db_ssp/final_scaffolds.functionalproteins.fasta | perl -pe "s/>/\n>/g" > $temp_ssp/$1_oneline.fasta
grep -A1 --no-group-separator -f $temp_ssp/$1-list1.txt $temp_ssp/$1_oneline.fasta > $ssp/$1-significant.fasta
echo "Number in fasta file"
grep -o ">" $ssp/$1-significant.fasta | wc -l
