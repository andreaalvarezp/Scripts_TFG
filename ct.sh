#!/bin/bash
#
#SBATCH --exclusive
#SBATCH --job-name=ct
#SBATCH --mem=200GB
#SBATCH -n 20
#SBATCH -p bigmem
#SBATCH --ntasks=1
#SBATCH --time=100:00:00

# Defino variables
resources="/home/aalvarez/bbmap/resources"
bbtool="/home/aalvarez/bbmap"
db="/home/aalvarez/db"
cap3="/home/aalvarez/CAP3"
ct="/home/aalvarez/Ct"
output="/home/aalvarez/Ct/output"
query="/home/aalvarez/Ct/SRA/query"
blast="/home/aalvarez/Ct/blast"
tables="/home/aalvarez/Ct/tables"
results3="/home/aalvarez/Ct/blast/$1"
results4="/home/aalvarez/Ct/blast/$1/NCBI"
results5="/home/aalvarez/Ct/blast/$1/aftercap"
BB="/home/aalvarez/Ct/BBTools"
dir3="/home/aalvarez/Ct/sambam"
dir=$PWD

module load BBMap/38.87-GCC-8.2.0-2.31.1
module load Trinity/2.11.0-foss-2019a-Python-3.7.2
module load DIAMOND/0.9.24-GCC-8.2.0-2.31.1
module load BWA/0.7.17-GCC-8.2.0-2.31.1
module load SAMtools/1.9-GCC-8.2.0-2.31.1

echo "Condición de muestras de Colletotrichum tofieldiae: $1"
echo "Posibles argumentos: invitro, plusP, minusP"
echo "Objetivo: Comprobación de la eficiencia del procesado en la identificación de micovirus en datos de RNAseq de distintos aislados de Coleltotrichum tofieldiae"

## 1. PRETRATAMIENTO CON BBTOOLS

# Sketch: me devuelve una tabla tabulada donde compruebo que las principales especies de mi muestra es Colletotrichum en el caso de las muestras in vitro y Arabidopsis en el caso de las muestras in planta
$bbtool/sendsketch.sh in=$query/Ct_$1.fastq.gz out=$BB/sketch_Ct_$1.txt reads=200000 merge printname0=f records=20 overwrite=true color=false depth depth2 volume sortbyvolume contam2=genus nt ow

# Clumpify: crea grupos y ordena las secuencias para facilitarsu manejo por las herramientas posteriores y reduce el tiempo de computación.
$bbtool/clumpify.sh -Xmx100G pigz=t unpigz=t zl=4 reorder in=$query/Ct_$1.fastq.gz out=$BB/Ct_$1_clump.fq.gz passes=1

# BBDuk:  elimina diversas zonas de las secuencias cogidas de librerías del programa utilizando unos K-mers de referencia de longitud escogida por el usuario.
# 1.1. Retiro adaptadores
$bbtool/bbduk.sh -Xmx75G -threads=15 ktrim=r ordered minlen=10 mink=10 tbo tpe rcomp=t overwrite=true k=15 hdist=1 hdist2=1 zl=4 ow=true in=$BB/Ct_$1_clump.fq.gz out=$BB/Ct_$1_temp2.fq.gz rqc=hashmap outduk=$output/Ct_$1_ktrim_kmerStats1.txt stats=$output/Ct_$1_ktrim_scaffoldStats1.txt loglog ref=$resources/adapters.fa

# 1.2.Retiro artefacts: ME SALTO ESTE PASO ME DA MAL
$bbtool/bbduk.sh -Xmx140G -threads=20 qtrim=r ordered k=25 hdist=1 hdist2=1 zl=4 cf=t barcodefilter=crash ow=true in=$BB/Ct_$1_temp2.fq.gz out=$BB/Ct_$1_temp3.fq.gz outm=$output/Ct_$1_Temp3_outmatch.fq.gz outduk=$output/Ct_$1_kmerStats2.txt stats=$output/Ct_$1_scaffoldStats2.txt loglog ref=$resources/phix_adapters.fa.gz,$resources/lambda.fa.gz,$resources/sequencing_artifacts.fa.gz

# 1.3. Retiro secuencias cortas
$bbtool/bbduk.sh -Xmx75G -Xms75G -threads=15 ordered overwrite=true k=20 hdist=1 zl=4 ow=true in=$BB/Ct_$1_temp2.fq.gz out=$BB/Ct_$1_temp4.fq.gz outm=$output/Ct_$1_Temp4_outmatch.fq.gz outduk=$output/Ct_$1_kmerStats3.txt stats=$output/Ct_$1_scaffoldStats3.txt loglog ref=$resources/short.fa

# 1.4. Retiro secuencias ribosómicas
$bbtool/bbduk.sh -Xmx75G -Xms75G -threads=15 ordered k=31 ref=$resources/ribokmers.fa.gz ow=true in=$BB/Ct_$1_temp4.fq out=$BB/Ct_$1_temp5.fq.gz

# Normalizacion de de la cobertura de los datos
$bbtool/bbnorm.sh -Xmx150G -Xms75G -threads=15 in=$BB/Ct_$1_temp5.fq.gz out=$BB/Ct_$1_clean.fq.gz

## 2. ENSAMBLAJE DE NOVO CON TRINITY

Trinity --seqType fq --max_memory 200G --CPU 40 --single $BB/Ct_$1_clean.fq.gz --output $ct/trinity_$1 --full_cleanup --verbose

## 3. PRIMER BLAST DIAMOND: BASE DE DATOS DE VIRUS CUSTOMIZADA

# 3.1. Creo base de datos a partir de las secuencias descargadas del NCBI
diamond makedb -p10 --in $db/virus_db.fasta --db $db/viraldb 

# 3.2. BLAST-X de las lecturas ensambladas contra la base de datos de virus personalizada. Se realiza dos veces para tener los resultados en dos formatos distintos.
diamond blastx -d $db/viraldb -q $ct/trinity_$1.Trinity.fasta -o $blast/align_trinity-$1-ViralaliNEW.blast -f 0 -k 2 -e 0.00001 --unal 0 --more-sensitive
diamond blastx -d $db/viraldb -q $ct/trinity_$1.Trinity.fasta -o $blast/align_trinity-$1-ViralNEW.blast -f 6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore -k 2 -e 0.00001 --unal 0 --more-sensitive

# Tabla de resultados del Blast after Trinity en Excel

# mkdir $results3/temp
grep -e ">" $blast/align_trinity-$1-ViralaliNEW.blast > $results3/list.txt
cut -d "[" -f 2 $results3/list.txt > $results3/list1.txt
sed -i 's/]/ /g' $results3/list1.txt
cut -d "[" -f 1 $results3/list.txt > $results3/list2.txt
sed -i 's/>//g' $results3/list2.txt
paste $blast/align_trinity-$1-ViralNEW.blast $results3/list2.txt > $results3/NCBItable1.blast
paste $results3/NCBItable1.blast $results3/list1.txt > $results3/NCBItable2.blast
cat <(head -1 header) $results3/NCBItable2.blast > $tables/$1_ViraltablefinalNEW.blast

# Extraigo los scaffolds correspondientes a cada elemento de la lista1

grep -e "Query=" $blast/align_trinity-$1-ViralaliNEW.blast > $results3/temp/Scaffoldslistnew.txt
cat $results3/temp/Scaffoldslistnew.txt | sort | uniq > $results3/temp/Scaffoldslistnewuni.txt
cut  -d " " -f 2 $results3/temp/Scaffoldslistnewuni.txt > $results3/$1-list1.txt
echo "Number of scaffolds in list1"
grep -o "TRINITY" $results3/$1-list1.txt | wc -l

# Extraigo archivo .fasta del output de Trinity 

awk '{if($0 ~ /^>/){print $0} else {printf $0}}' $ct/trinity_$1.Trinity.fasta | perl -pe "s/>/\n>/g" > $ct/trinity_$1_oneline.fasta
grep -A1 --no-group-separator -f $results3/$1-list1.txt $ct/trinity_$1_oneline.fasta > $ct/$1-significant1NEW.fasta
echo "Numero de secuencias en archivo fasta"
grep -o ">" $ct/$1-significant1NEW.fasta | wc -l

# A partir de aquí se siguen dos procesados distintos para posteriormente comparar los resultados

### 4. PROCESADO 1

## 4.1. BLAST VS NCBI CON OUTPUT DE TRINITY

# 4.1.1. Creo una base de datos con las secuencias del NCBI descargadas, la NCBI nr (no-redundat). Dicha base de datos fue descargada en marzo de 2021.
diamond makedb --in nr.gz -d $db/NEWNR/nr

# 4.1.2. BLAST-X de las lecturas ensambladas contra la base de datos del NCBI. Se realiza dos veces para tener los resultados en dos formatos distintos.
echo "Blast contra NCBI Ct_$1"
diamond blastx -p20 -d $db/NEWNR/nr -q $ct/$1-significant1NEW.fasta -o $blast/$1-NCBI_afterTrinity.blast -f 0 -k 1 -e 0.00001 --unal 0 --more-sensitive
diamond blastx -d $db/NEWNR/nr -q $ct/$1-significant1NEW.fasta -o $blast/$1-NCBItab_afterTrinity.blast -f 6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore -k 2 -e 0.001 --unal 0 --more-sensitive

# Tabla de resultados del blast contra NCBI a Excel

grep -e ">" $blast/$1-NCBI_afterTrinity.blast > $results4/list_af.txt
cut -d "[" -f 2 $results4/list_af.txt > $results4/list1_af.txt
sed -i 's/]/ /g' $results4/list1_af.txt
cut -d "[" -f 1 $results4/list_af.txt > $results4/list2_af.txt
sed -i 's/>//g' $results4/list2_af.txt
paste $blast/$1-NCBItab_afterTrinity.blast $results4/list2_af.txt > $results4/$1-NCBItable1_af.blast
paste $results4/$1-NCBItable1_af.blast $results4/list1_af.txt > $results4/$1-NCBItable2_af.blast
cat <(head -1 header) $results4/$1-NCBItable2_af.blast > $tables/$1-NCBItablefinal_af.blast

### 5. PROCESADO 2

# 5.1. REENSAMBLAJE CON CAP3

echo "CAP3 processing $1"
mkdir $dir/cap3_data_$1
cp $ct/$1-significant1NEW.fasta $ct/cap3_data/$1/
$cap3/cap3 $ct/cap3_data/$1/$1-significant1NEW.fasta > $ct/cap3_data/$1/$1-significant1NEW.cap3
cat $dir/cap3_data/$1/*cap.singlets $dir/cap3_data/$1/*cap.contigs > $ct/$1_mapNEW.fasta

# 5.2. SEGUNDO BLAST DIAMOND TRAS CAP3 contra la base de datos de virus customizada como segundo filtro de los datos. También se realiza dos veces para dar lugar a la salida en dos formatos distintos

echo "Blast contra base de virus customizada"
diamond blastx -p20 -d $db/viraldb -q $ct/$1_mapNEW.fasta -o $blast/$1-ViralaliNEWafterCAP.blast -f 0 -k 1 -e 0.00001 --unal 0 --more-sensitive
diamond blastx -p20 -d $db/viraldb -q $ct/$1_mapNEW.fasta -o $blast/$1-ViralNEWafterCAP.blast -f 6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore -k 1 -e 0.00001 --unal 0 --more-sensitive

# Tabla de resultados del Blast after CAP3 en Excel

grep -e ">" $blast/$1-ViralaliNEWafterCAP.blast > $results4/list.txt
cut -d "[" -f 2 $results4/list.txt > $results4/list1.txt
sed -i 's/]/ /g' $results4/list1.txt
cut -d "[" -f 1 $results4/list.txt > $results4/list2.txt
sed -i 's/>//g' $results4/list2.txt
paste $blast/$1-ViralNEWafterCAP.blast $results4/list2.txt > $results4/NCBItable1.blast
paste $results4/NCBItable1.blast $results4/list1.txt > $results4/NCBItable2.blast
cat <(head -1 header) $results4/NCBItable2.blast > $tables/$1_ViraltablefinalNEWafterCAP.blast

# Extraigo los scaffolds correspondientes a cada elemento de la lista2

grep -e "Query=" $blast/$1-ViralaliNEWafterCAP.blast > $results4/$1-list2.txt
cat $results4/$1-list2.txt | sort | uniq > $results4/$1-list2uni.txt
cut -d " " -f 2 $results4/$1-list2uni.txt > $results4/$1-list2.txt
echo "Numero de scaffolds en la list2"
grep -o -e "TRINITY" -e "Contig" $results4/$1-list2.txt | wc -l

# Extraigo archivos fasta de los resultados de CAP3

awk '{if($0 ~ /^>/){print $0} else {printf $0}}' $ct/$1_mapNEW.fasta | perl -pe "s/>/\n>/g" > $ct/$1_mapNEW_oneline.fasta
grep -A1 --no-group-separator -f $results4/$1-list2.txt $ct/$1_mapNEW_oneline.fasta > $ct/$1-significant2NEW.fasta
echo "Numero de secuencias en archivo $1-significant2NEW.fasta"
grep -o ">" $ct/$1-significant2NEW.fasta | wc -l 

## 5.3. BLAST TO NCBI CON OUTPUT DE CAP3: Identifico las secuencias virales con este blastx

echo "BLAST to NCBI con output de CAP3 $1"
diamond blastx -p20 -d $db/NEWNR/nr -q $ct/$1-significant2NEW.fasta -o $blast/$1-NCBI.blast -f 0 -k 1 -e 0.00001 --unal 0 --more-sensitive
diamond blastx -d $db/NEWNR/nr -q $ct/$1-significant2NEW.fasta -o $blast/$1-NCBItab.blast -f 6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore -k 1 -e 0.00001 --unal 0 --more-sensitive

# Tabla de resultados del segundo blast a Excel

grep -e ">" $blast/$1-NCBI.blast > $results5/list.txt
cut -d "[" -f 2 $results5/list.txt > $results5/list1.txt
sed -i 's/]/ /g' $results5/list1.txt
cut -d "[" -f 1 $results5/list.txt > $results5/list2.txt
sed -i 's/>//g' $results5/list2.txt
paste $blast/$1-NCBItab.blast $results5/list2.txt > $results5/$1-NCBItable1.blast
paste $results5/$1-NCBItable1.blast $results5/list1.txt > $results5/$1-NCBItable2.blast
cat <(head -1 header) $results5/$1-NCBItable2.blast > $tables/$1-NCBItablefinal.blast
