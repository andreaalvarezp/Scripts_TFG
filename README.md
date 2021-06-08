# README: Búsqueda de motivos genómicos con posible relevancia biológica en tres aislados de hongos del género *Plectosphaerella*

**Grado en Biotecnología, Universidad Politécnica de Madrid**

**Autor: Andrea Álvarez Pérez**

**Curso 2020-2021**

Este Github recoge los scripts bioinformáticos utilizados para la realización del **Trabajo de Fin de Grado de Biotecnología** de la alumna Andrea Álvarez Pérez. Cada uno de ellos está destinado a una tarea concreta y todos los scripts, tanto en R como en bash, están comentados con la función y el objetivo de cada una de las órdenes y módulos ejecutados.

De igual manera, en este archivo se explica de manera concisa el uso de cada uno de los script y se clasifican según la sección del Trabajo de Fin de Grado a la que pertenecen. El uso y los parámetros de cada módulo se especifican dentro de cada uno de los script.

**Importante**: si se quiere hacer uso de los scripts es obligatorio modificar los *paths* establecidos en el primer bloque de cada script y redefinir las variables de manera personalizada para cada usuario.

## 1. BÚSQUEDA Y DETERMINACIÓN DE MICOVIRUS EN *PLECTOSPHAERELLA*

### 1.1. *Plec.sh*

Este módulo aplica un procesa bioinformático a los datos de RNA-seq brutos de *Plectosphaerella* que permite la búsqueda e identificación de micovirus en datos de RNAseq de distintos aislados de *Plectospherella*. Uso:

```bash
$ sbatch Plec.sh $1
```

siendo ``$1`` la cepa de *Plectosphaerella* a analizar. Posibles argumentos:
- PcBMM
- Pc2127
- P0831

Este script aplica diversos módulos cargados desde el cluster de supercomputación del CBGP:

```bash
module load BBMap/38.87-GCC-8.2.0-2.31.1
module load Trinity/2.11.0-foss-2019a-Python-3.7.2
module load DIAMOND/0.9.24-GCC-8.2.0-2.31.1
```

### 1.2. *ct.sh*

Este módulo es idéntico al anterior, pero posee las variables redefinidas para utilizarlas con datos del hongo *Colletotrichum tofieldiae*. Sirve para la búsqueda e identificación de micovirus en datos de RNAseq. Uso:

```bash
$ sbatch ct.sh $1
```

siendo ``$1`` las condiciones de las muestras de *Colletotrichum tofieldiae*. Posibles argumentos:
- *invitro*
- plusP
- minusP

Los módulos utilizados fueron los mismos que con *Plectosphaerella*.

### 1.3. *gilbert.sh*

Este módulo aplica el procesa bioinformático llevado a cabo por Gilbert, K., *et al*, 2019. Se aplica para comprobar la eficiencia del procesado llevado a cabo anteriormente y comparar los resultados obtenidos. Uso:

```bash
$ sbatch gilbert.sh $1
```

siendo ``$1`` las condiciones de las muestras de *Colletotrichum tofieldiae*. Posibles argumentos:
- *invitro*
- plusP
- minusP

Este script aplica diversos módulos cargados desde el cluster de supercomputación del CBGP:

```bash
module load BBMap/38.87-GCC-8.2.0-2.31.1
module load Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1
module load Trinity/2.11.0-foss-2019a-Python-3.7.2
module load DIAMOND/0.9.24-GCC-8.2.0-2.31.1
module load HMMER/3.2.1-GCC-8.2.0-2.31.1
module load Miniconda3/4.7.10
```

Además requiere de la creación y activación de un entorno virtual para ejecutar **cutadapt** y **TransDecoder**:

```bash
conda create -n AndreaConda python=3.7 anaconda
source activate AndreaConda
conda install -c bioconda cutadapt
# conda install -c bioconda transdecoder
```

## 2. CANDIDATOS DE SECUENCIAS DE *PLECTOSPHAERELLA* COMO MIMÉTICAS DE DOS FAMILIAS DE PÉPTIDOS FITOREGULADORES DE *ARABIDOPSIS*

### 2.1. *scoop.sh*

Este módulo se utilizó para la búsqueda de motivos SCOOP y SSP en el proteoma de *Plectosphaerella*. Uso:

```bash
$ sbatch scoop.sh $1
```

siendo ``$1`` la secuencia que se quiere utilizar como *query* para buscar en los proteomas del hongo. Posibles argumentos:
- scoop10
- scoop12
- fusarium
- verticilium
- magnaporthe
- proscoop

En cada uno de los archivos con el nombre del argumento deben estar alojadas la/s secuencia/s correspondientes a cada categoría. El script aplica un BLAST-P para su búsqueda y luego obtiene a partir de los resultados el archivo FASTA correspondiente que será utilizado como entrada en MEME. El módulos utilizado es:

```bash
module load BLAST+/2.9.0-gompi-2019a
```

### 2.2. *SSP_heatmaps.R*

Módulo utilizado para la obtención de mapas de calor de expresión de genes de *Plectosphaerella* interesantes en la identificación de miméticos de péptidos SSP. Es el único de los scripts que se ha programado con R, haciendo uso de diversas librerías:

```R
library(gplots)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
```

Cada una de estas librerías se utilizó con un uso concreto:
- ``gplots`` se utilizó para la representación de diagramas y dendogramas. 
- ``dplyr`` sirvió para el cálculo de las matrices de distancias y la agrupación de los datos en clústeres.
- ``ComplexHeatMap`` permitió la representación de los mapas de calor. Se hizo uso de múltiples de sus opciones y parámetros, comentados en el código.
- ``RColorBrewer`` permitió crear un código de color para cada cluster.

## 3. BÚSQUEDA Y DETERMINACIÓN DE FAMILIAS DE CAZymas EN *PLECTOSPHAERELLA*

### 3.1. *CAZymas.sh*

Este módulo ejecuta múltiples tareas:
1. Eejecuta dbCAN2 para la determinación de genes ue codifican para CAZymas en los proteomas y genomas de *Plectosphaerella*.
3. Divide los resultados según la familia o grupo al que pertenecen (AA, CE, CBM, GH, GT, PL).
4. Prepara el archivo FASTA correspondiente para su posterior ejecución en SECRETOOL.
5. Alinea las secuencias de las tres cepas del hongo que pertenecen a una familia con ClustalOmega.
6. Genera un árbol filogenético con IQTREE.

```bash
$ sbatch CAZymas.sh $1 $2
```

siendo ``$1`` la cepa de *Plectosphaerella* a analizar. Posibles argumentos:
- PcBMM
- Pc2127
- P0831

siendo ``$2`` la familia de CAZymas que queremos posteriormente filtrar. Posibles argumentos:
- AA
- CE
- CBM (no se considera familia, agrupa genes que contienen dominios de unión a carbohidratos)
- GH
- GT
- PL

Para la ejecución de este script es necesario crear un entorno virtual para poder utilizar dbCAN2 de manera remota:

```bash
module load Miniconda3/4.7.10
conda create -n run_dbcan python=3.8 diamond hmmer prodigal -c conda-forge -c bioconda
source activate run_dbcan
pip install run-dbcan==2.0.11
conda install -c bioconda run-dbcan==2.0.11
```

Además de la ejecución de diversos módulos:

```bash
module load DIAMOND/0.9.24-GCC-8.2.0-2.31.1
module load BLAST+/2.9.0-gompi-2019a
module load Clustal-Omega/1.2.4-GCC-8.2.0-2.31.1
module load IQ-TREE/1.6.12-foss-2019a
```

Este script da como resultado seis archivos, uno por cada grupo de CAZymas, con su correpondiente FASTA, alineamiento y archivo *.tree* con el árbol que será posteriomente representado con iTOL.
