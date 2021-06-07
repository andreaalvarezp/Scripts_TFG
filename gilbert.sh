#!/bin/bash
#
#SBATCH --job-name=ct
#SBATCH --mem=220GB
#SBATCH -n 20
#SBATCH -p bigmem
#SBATCH --ntasks=1
#SBATCH --time=200:00:00

# Defino variables
resources="/home/aalvarez/bbmap/resources"
bbtool="/home/aalvarez/bbmap"
db="/home/aalvarez/db"
ct="/home/aalvarez/Ct_gilbert"
BB="/home/aalvarez/Ct_gilbert/BBTools"
output="/home/aalvarez/Ct_gilbert/output"
query="/home/aalvarez/Ct_gilbert/paired-query"
blast="/home/aalvarez/Ct_gilbert/blast"
tables="/home/aalvarez/Ct_gilbert/tables"
bowtie2="/home/aalvarez/Ct_gilbert/bowtie2"
cutadapt="/home/aalvarez/Ct/cutadapt"
trinity="/home/aalvarez/Ct_gilbert/trinity"
transdecoder="/home/aalvarez/Ct_gilbert/transdecoder"
HMMscan="/home/aalvarez/Ct_gilbert/HMMscan"
dir=$PWD

## Cargo módulos y activo el entorno virtual

module load BBMap/38.87-GCC-8.2.0-2.31.1
module load Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1
module load Trinity/2.11.0-foss-2019a-Python-3.7.2
module load DIAMOND/0.9.24-GCC-8.2.0-2.31.1
module load Miniconda3/4.7.10
conda create -n AndreaConda python=3.7 anaconda
source activate AndreaConda
conda install -c bioconda cutadapt

echo "Artículo Gilbert pipeline: $1"
echo "Posibles argumentos: invitro, minusP, plusP"
echo "Objetivo: comprobar la fiabilidad del pipeline de procesado propuesto repitiendo el pipeline llevado a cabo por Gilbert, K., et al. (2019)"

## 1. SEPARAR PAIRED-ENDS EN DOS ARCHIVOS DISTINTOS

$bbtool/reformat.sh in=$ct/SRA/Ct_$1_$2.fastq.gz out1=$ct/SRA/$1_$2_R1.fq.gz out2=$ct/SRA/$1_$2_R2.fq.gz
$bbtool/reformat.sh in=$ct/SRA/$1Ct_$2.fastq.gz out1=$ct/SRA/$1_$2_R1.fq.gz out2=$ct/SRA/$1_$2_R2.fq.gz

## 2. QUITO ADAPTADORES CON CUTADAPT

module load Miniconda3/4.7.10
source activate AndreaConda

cutadapt -o $cutadapt/$1_1.fq.gz -p $cutadapt/$1_2.fq.gz $paired/$1_R1.fq.gz $paired/$1_R2.fq.gz

# Retirar los adaptadores con cutadapt sin conocer los parámetros exactos del trabajo de Gilbert dio lugar a archivos del mismo tamaño.
# ALTERNATIVA: QUITARLOS CON BBTOOLS. Solo se llevaron a cabo los pasos de agrupación para reducir el tamaño de los archivos y el tiempo de computación posterior, así como la eliminación de los adapatadores con BBDuk.

# Sketch
$bbtool/sendsketch.sh in=$query/Ct_$1.fastq.gz out=$BB/sketch_Ct_$1.txt reads=200000 merge printname0=f records=20 overwrite=true color=false depth depth2 volume sortbyvolume contam2=genus nt ow

# Clumpify
$bbtool/clumpify.sh -Xmx100G pigz=t unpigz=t zl=4 reorder in1=$paired/$1_R1.fq.gz in2=$paired/$1_R2.fq.gz out1=$BB/$1_clump1.fq.gz out2=$BB/$1_clump2.fq.gz passes=1

# BBDuk
# 1.1. Retiro adaptadores
$bbtool/bbduk.sh -Xmx75G -threads=15 ktrim=r ordered minlen=10 mink=10 tbo tpe rcomp=t overwrite=true k=15 hdist=1 hdist2=1 zl=4 ow=true in1=$BB/$1_clump1.fq.gz in2=$BB/$1_clump2.fq.gz out1=$BB/$1_temp1.fq.gz out2=$BB/$1_temp2.fq.gz rqc=hashmap outduk=$BB/output/$1_ktrim_kmerStats.txt stats=$BB/output/$1_ktrim_scaffoldStats.txt loglog ref=$resources/adapters.fa

## 3. BOWTIE2

# 3.1. CREAR INDEX DEL GENOMA DE CT: construye un índice Bowtie a partir de un conjunto de secuencias de DNA, son todo lo que se necesita para alinear las lecturas con esa referencia.

bowtie2-build -f $genome/Colletotrichum_genome.fna $bowtie2/Ct_bt2base # El genoma de Colletotrichum se obtuvo del NCBI

# 3.2. BOWTIE2 - ALINEAMIENTO: Se utilizó la opción --un para poder conservar aquellas secuencias que no alineaban.

bowtie2 -x $bowtie2/genome/Ct_bt2base -1 $BB/$1_temp1.fq.gz -2 $BB/$1_temp2.fq.gz --un-conc-gz $bowtie2/$1_align.fq.gz

# 4. ENSAMBLAJE DE NOVO CON TRINITY

echo "Sample $1: ensamblaje de novo con Trinity"

Trinity --seqType fq --max_memory 200G --CPU 40 --left $bowtie2/$1_align.fq.1.gz --right $bowtie2/$1_align.fq.2.gz --output $ct/trinity_out2 --full_cleanup --verbose
