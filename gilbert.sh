#!/bin/bash
#
#SBATCH --job-name=gilbert
#SBATCH --mem=18GB
#SBATCH -n 20
#SBATCH -p medium
#SBATCH --ntasks=1
#SBATCH --time=20:00:00

# Defino variables
resources="/home/aalvarez/bbmap/resources"
bbtool="/home/aalvarez/bbmap"
db="/home/aalvarez/db"
ct="/home/aalvarez/Ct_gilbert"
query="/home/aalvarez/Ct_gilbert/paired-query"
output="/home/aalvarez/Ct_gilbert/output"
blast="/home/aalvarez/Ct_gilbert/blast"
tables="/home/aalvarez/Ct_gilbert/tables"
bowtie2="/home/aalvarez/Ct_gilbert/bowtie2"
cutadapt="/home/aalvarez/Ct/cutadapt"
transdecoder="/home/aalvarez/Ct_gilbert/trinity_$1.Trinity.fasta.transdecoder_dir"
hmmer="/home/aalvarez/Ct_gilbert/hmmer"
dir=$PWD

## Cargo módulos y activo el entorno virtual

module load BBMap/38.87-GCC-8.2.0-2.31.1
module load Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1
module load Trinity/2.11.0-foss-2019a-Python-3.7.2
module load DIAMOND/0.9.24-GCC-8.2.0-2.31.1
module load HMMER/3.2.1-GCC-8.2.0-2.31.1
module load Miniconda3/4.7.10
conda create -n AndreaConda python=3.7 anaconda
source activate AndreaConda
conda install -c bioconda cutadapt
conda install -c bioconda transdecoder

echo "Artículo Gilbert pipeline: $1"
echo "Posibles argumentos: invitro, minusP, plusP"
echo "Objetivo: comprobar la fiabilidad del pipeline de procesado propuesto repitiendo el pipeline llevado a cabo por Gilbert, K., et al. (2019)"

## 1. SEPARAR PAIRED-ENDS EN DOS ARCHIVOS DISTINTOS

$bbtool/reformat.sh in=$ct/SRA/Ct_$1_$2.fastq.gz out1=$ct/SRA/$1_$2_R1.fq.gz out2=$ct/SRA/$1_$2_R2.fq.gz

## 2. QUITO ADAPTADORES CON CUTADAPT

cutadapt -o $cutadapt/$1_1.fq.gz -p $cutadapt/$1_2.fq.gz $query/$1_R1.fq.gz $query/$1_R2.fq.gz

# Retirar los adaptadores con cutadapt sin conocer los parámetros exactos del trabajo de Gilbert dio lugar a archivos del mismo tamaño.

## 3. BOWTIE2

# 3.1. CREAR INDEX DEL GENOMA DE CT: construye un índice Bowtie a partir de un conjunto de secuencias de DNA, son todo lo que se necesita para alinear las lecturas con esa referencia.

bowtie2-build -f $genome/Colletotrichum_genome.fna $bowtie2/Ct_bt2base # El genoma de Colletotrichum se obtuvo del NCBI

# 3.2. BOWTIE2 - ALINEAMIENTO: Se utilizó la opción --un para poder conservar aquellas secuencias que no alineaban.

bowtie2 -x $bowtie2/genome/Ct_bt2base -1 $query/$1_R1.fq.gz -2 $query/$1_R2.fq.gz --un-conc-gz $bowtie2/$1_align.fq.gz

# MODIFICACIÓN PIPELINE: NORMALIZACIÓN

$bbtool/bbnorm.sh -Xmx150G -Xms75G -threads=15 in1=$bowtie2/$1_align.1.fq.gz in2=$bowtie2/$1_align.2.fq.gz out1=$bowtie2/$1_align_norm.1.fq.gz out2=$bowtie2/$1_align_norm.2.fq.gz

# 4. ENSAMBLAJE DE NOVO CON TRINITY

echo "Sample $1: ensamblaje de novo con Trinity"

Trinity --seqType fq --max_memory 200G --CPU 40 --left $bowtie2/$1_align.fq.1.gz --right $bowtie2/$1_align.fq.2.gz --output $ct/trinity_$1 --full_cleanup --verbose

# 6. IDENTIFICACIÓN DE REGIONES CODIFICANTES CON TRANSDECODER

# Para la instalación de TransDecoder es necesario hacer uso del entorno virtual

TransDecoder.LongOrfs -t trinity_$1.Trinity.fasta

# 7. IDENTIFICACIÓN DE SECUENCIAS CON SIMILITUD LIMITADA A DOMINIOS DE PROTEÍNAS VIRALES CON HMMER

# 7.1. Descarga de la base de datos de Pfam

wget -O $ct/db/Pfam-A.hmm.gz "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"

# 7.2. Busco los perfiles HMM que me interesan: "RdRP_1, RdRP_2, RdRP_3, RdRP_4, RdRP_5, Mitovir_RNA_pol"

hmmfetch $ct/db/Pfam-A.hmm.gz $1 > $1.hmms

# Junto todos los archivos en uno solo

cat RdRP_*.hmms Mitovir_RNA_pol.hmms > RdRP.hmm

# 7.3. Creo la base de datos con hmmpress

hmmpress RdRP.hmm

# 7.4. Busco las proteínas con hmmscan

hmmscan RdRP.hmm $transdecoder/longest_orfs.pep > $hmmer/hmmer_$1

# 8. COMPROBACIÓN DE LOS RESULTADOS CON BLAST
# 8.1. Creo un archivo oneline de los ensamblajes de Trandecoder y creo una lista manual con los identificadores de los hit de HMMER
# Con ese resultado obtengo el archivo FASTA con las proteínas que han dado hit

awk '{if($0 ~ /^>/){print $0} else {printf $0}}' $transdecoder/longest_orfs.pep | perl -pe "s/>/\n>/g" > $transdecoder/longest_orfs_oneline.fasta
grep -A1 --no-group-separator -f $hmmer/list1_$1.txt $transdecoder/longest_orfs_oneline.fasta > $hmmer/hmmer-$1.fasta

# 9. COMPROBACIÓN DE RESULTADOS: BLAST CONTRA NCBI nr

diamond blastp -d $db/NEWNR/nr -q $ct/hmmer/hmmer-$1.fasta -o $ct/hmmer-NCBI_$1.blast -f 0 -k 2 -e 0.001 --unal 0 --more-sensitive
