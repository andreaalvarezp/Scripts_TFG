#!/bin/bash
#
#SBATCH --job-name=dbCAN2
#SBATCH --mem=100GB
#SBATCH -n 20
#SBATCH -p test
#SBATCH --ntasks=1
#SBATCH --time=10:00

## Defino variables

caz="/home/aalvarez/CAZymas"
genoma="/home/aalvarez/CAZymas/genomas"
output="/home/aalvarez/CAZymas/output_dbCAN"
proteoma="/home/aalvarez/CAZymas/proteomas"
tmp="/home/aalvarez/CAZymas/tmp"
fasta="/home/aalvarez/CAZymas/fasta"
tmp2="/home/aalvarez/CAZymas/tmp_genoma"
fasta2="/home/aalvarez/CAZymas/fasta_genoma"
secretool="/home/aalvarez/CAZymas/SECRETOOL_proteoma"
cupp="/home/aalvarez/CAZymas/CUPP-dbCAN"
prot="/home/aalvarez/CAZymas/output_dbCAN/proteoma_filtered"
final="/home/aalvarez/CAZymas/output_dbCAN/proteoma_filtered/fasta-final"
fam="/home/aalvarez/CAZymas/output_dbCAN/proteoma_filtered/fasta-final/$1"

## Cargo módulos

module load DIAMOND/0.9.24-GCC-8.2.0-2.31.1
module load BLAST+/2.9.0-gompi-2019a

echo "dbCAN2: automated CAZyme annotation cepa $1 Plectosphaerella"
echo "Posibles argumentos: PcBMM, Pc2127, P0831"
echo "Objetivo: identificación de genes en el proteoma y el genoma de Plectosphaerella que codifican para CAZymas"

## 1. INSTALACIÓN DE dbCAN2

# 1.1. Creo entorno virtual
module load Miniconda3/4.7.10
conda create -n run_dbcan python=3.8 diamond hmmer prodigal -c conda-forge -c bioconda
source activate run_dbcan
pip install run-dbcan==2.0.11
conda install -c bioconda run-dbcan==2.0.11

# 1.2. Descarga de base de datos
test -d db || mkdir db
cd db \
    && wget http://bcb.unl.edu/dbCAN2/download/CAZyDB.07312019.fa.nr && diamond makedb --in CAZyDB.07312019.fa.nr -d CAZy \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V8.txt && mv dbCAN-HMMdb-V8.txt dbCAN.txt && hmmpress dbCAN.txt \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-1.hmm && hmmpress tf-1.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-2.hmm && hmmpress tf-2.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/stp.hmm && hmmpress stp.hmm \
    && cd ../ && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.fna \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.faa \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.gff

# 1.3. Instalación de SignalP
mkdir -p run_dbcan/tools && 
cd run_dbcan/tools/ && tar xzf signalp-4.1g.Linux.tar.gz && cd signalp-4.1

# 1.4. Edito el archivo signalp con los directorios del ejecutable y el directorio de salida

# 1.5. Instalación
cd run_dbcan/tools/signalp-4.1
sudo cp signalp /usr/bin/signalp
sudo chmod 755 /usr/bin/signalp

# 1.6. Compruebo el programa
run_dbcan.py EscheriaColiK12MG1655.fna prok --out_dir output_EscheriaColiK12MG1655

## 2. dbCAN2 PROTEOMA PLECTOSPHAERELLA

run_dbcan.py $proteoma/$1_functional_proteins.fasta protein --out_dir $output/dbCAN_$1_proteoma
echo "dbCAN realizado con exito"

# 2.1. Filtro manual, me quedo con los que detectan al menos dos de las herramientas

# 2.2. Obtengo archivo FASTA

grep -e "CBGP_" $output/dbCAN_$1_proteoma/overview.txt > $tmp/list_$1.txt
cut -f 1 $tmp/list_$1.txt > $tmp/$1_proteinID.txt
awk '{if($0 ~ /^>/){print $0} else {printf $0}}' $proteoma/$1_functional_proteins.fasta | perl -pe "s/>/\n>/g" > $tmp/$1_oneline.fasta
grep -A1 --no-group-separator -f $tmp/$1_proteinID.txt $tmp/$1_oneline.fasta > $fasta/dbCAN_$1.fasta

## 3. dbCAN2 GENOMA PLECTOSPHAERELLA

run_dbcan.py $genoma/$1_genome_assembly.fasta meta --out_dir $output/dbCAN_$1_genoma
echo "dbCAN realizado con exito"

# 3.2. Obtengo archivo FASTA

grep -e "Scaffold" $output/dbCAN_$1_genoma/overview.txt > $tmp2/list_$1.txt
cut -d "_" -f 1,2,3,4,5 $tmp2/list_$1.txt > $1_proteinID.txt
awk '{if($0 ~ /^>/){print $0} else {printf $0}}' $genoma/$1_genome.fasta | perl -pe "s/>/\n>/g" > $tmp2/$1_oneline_genoma.fasta
cat $tmp2/$1_proteinID.txt | sort | uniq > $tmp2/$1_proteinID_uniq.txt
grep -A1 --no-group-separator -f $tmp2/$1_proteinID_uniq.txt $tmp2/$1_oneline_genoma.fasta > $fasta2/dbCAN_$1.fasta

### MODIFICACIÓN DE LOS ARCHIVOS DE SALIDA

echo "Separación por familias y preparación de FASTA de output de dbCAN2"
echo "Los archivos sobre los que se trabaja se han filtrado manualmente para quedarnos con las CAZymas detectadas por al menos dos de las herramientas de dbCAN2"

echo "Aislado $1 - Familia $2"
echo "Posibles primeros argumentos: PcBMM, Pc2127, P0831"
echo "Posibles segundos argumentos: AA, CE, CBM, GH, GT, PL"

# 1. Me quedo con las líneas que tenga la familia $2
grep -n "$2" $prot/$1_filtered.txt > $prot/$1/$1_$2.txt

# 2. Me quedo solo con los identificadores de las proteinas
cut -f 1 $prot/$1/$1_$2.txt | cut -d ":" -f 2 > $prot/$1/$2_ID.txt

# 3. Identificadores únicos
cat $prot/$1/$2_ID.txt | sort | uniq > $prot/$1/$2_uniqID.txt
echo "Numero de CAZymas de la familia $2 en el aislado $1:" 
cat $prot/$1/$2_uniqID.txt | wc -l

# 4. Genero archivo oneline del proteoma
awk '{if($0 ~ /^>/){print $0} else {printf $0}}' $proteoma/$1_functional_proteins.fasta | perl -pe "s/>/\n>/g" > $prot/$1_oneline.fasta

# 5. Genero archivo FASTA
grep -A1 --no-group-separator -f $prot/$1/$2_uniqID.txt $prot/$1_oneline.fasta > $prot/fasta-final/$2/$1_$2.fasta

echo "Archivo FASTA generado con exito"

# 6. Modifico los identificadores sustituyendo el encabezado "CBGP" por el nombre de la cepa:
sed -i 's/CBGP/$1/g' "$prot/fasta-final/$2/$1_$2.fasta"

# 7. Sustituyo la coletilla "RA" por el nombre de la familia de CAZymas
sed -i 's/-RA prot/-$2 prot/g' "$prot/fasta-final/$2/$1_$2.fasta"

## A continuación hay unos pasos de modificación manual de los archivos de manera que debemos tener un archivo por cada familia de CAZymas con las secuencias de las tres cepas
## Es decir, manualmente se deben juntar los tres archivos "$1_$2.fasta" de cada familia, teniendo un total de 6 archivos (AA, CE, CBM, GH, GT, PL)
## Archivos necesarios: AA.fasta, CE.fasta, CBM.fasta, GH.fasta, GT.fasta y PL.fasta

module load Clustal-Omega/1.2.4-GCC-8.2.0-2.31.1
module load IQ-TREE/1.6.12-foss-2019a

echo "Construccion arbol filogenetico familia $2"

## 1. ALINEAMIENTO: utilizamos la última version de Clustal Omega para realizar el alineamiento múltiple de las secuencias

clustalo -i $fam/$2.fasta -o $fam/$2_alin.fasta -v
clustalo --full --percent-id -i $fam/$2_alin.fasta -v -o $fam/$2_matrix.csv

## 2. CONSTRUCCIÓN DEL ÁRBOL FILOGENÉTICO
# Metodo de sustitucion TEST: el programa hace un TEST para determinar cual es el mejor modelo de sustitución que se ajusta a mis datos y construye el árbol aplicándolo. 
# Para aportar robustez a lo árboles se estableció un valor de bootstrap de 1000.

iqtree -s $fam/$2_alin.fasta -nt AUTO -bb 1000 -m TEST
