# README: Búsqueda de motivos genómicos con posible relevancia biológica en tres aislados de hongos del género *Plectosphaerella*

**Grado en Biotecnología, Universidad Politécnica de Madrid**

**Autor: Andrea Álvarez Pérez**

**Tutores: Soledad Sacristán Benayas, Julio Luis Rodríguez Romero**

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

El esquema del procesado seguido se muestra en el link:

![ScreenShot](https://github.com/andreaalvarezp/Scripts_TFG/tree/main/Pipelines/pipeline_Plec.png)

A continuación se realiza una breve descripción de las utilidades que le damos a cada uno de los módulos:

#### BBTools (Bushnell, B., *et al,* 2017) [22]

Herramienta bioinformática multiproceso de la que vamos a ejecutar diversos módulos:

**1) Sendsketch** hace un boceto de mi ensamblaje, que se limitó a los primeros 200.000 reads. Compara los croquis de consulta creados a partir de los archivos FASTQ de entrada con los croquis de referencia alojados en un servidor remoto a través de Internet.

**2) Clumpify** es el módulo encargado de agrupar rápidamente las lecturas superpuestas en grupos. Se utiliza como una forma de aumentar la compresión de archivos y para acelerar el ensamblaje.

**3)	BBDuk** se utiliza para retirar adaptadores de las secuencias de Illumina; artefacts, que son repeticiones de aminoácidos o secuencias que no son necesarias; secuencias cortas y ribosómicas. 

**4)	BBNorm** está diseñado como herramienta de normalización de la cobertura mediante el muestreo descendente de las lecturas para dar una distribución de cobertura plana. Se busca minimizar el ruido técnico introducido en los datos durante el proceso de secuenciación con el fin de volverlos comparables entre sí, lo que acelera y mejora la calidad del ensamblaje posterior.

#### Trinity (Grabherr, M.G., *et al,* 2011) [23]

Trinity es un software que combina tres módulos independientes que se ejecutan de manera secuencial para procesar grandes volúmenes de lecturas de RNAseq.

Primero reconstruye los contigs a partir del set de lecturas, agrupa los superpuestos y los conectan para construir gráficos de Bruijin para cada componente. Un gráfico de Bruijin está definido por nodos donde cada uno de ellos representa una secuencia de una longitud fija de k nucleótidos (k-mer), donde k es considerablemente más corto que la longitud de la lectura (Grabherr, M.G., et al 2011), y permiten enumerar todas las posibles soluciones por las se pueden reconstruir secuencias lineales, por lo que cada ruta en el gráfico representa una posible transcripción. Por último, reconstruye los transcritos lineales resolviendo los gráficos de Bruijin individuales con las lecturas originales y puntúa las rutas basándose en procesos de programación dinámica.

#### DIAMOND (Buchfink, B., *et al,* 2015) [24]

DIAMOND es un software alineamiento de secuencias para búsquedas de proteínas y DNA traducido, diseñado para el análisis de alto rendimiento de datos de grandes secuencias. Es un programa de código abierto cuatro veces más rápido que BLAST-X en la comparación de lecturas cortas de DNA con la base de datos NCBI nr (non-redundant) manteniendo un nivel comparable de sensibilidad en alineaciones con un valor E (evalue) < 10-3.

#### CAP3 - Contig Assembly Program (Huang, X., 1999) [25]

Programa utilizado para el reensamblaje de contigs. Se trata de una herramienta que utiliza un filtro para eliminar pares de fragmentos que posiblemente no podrían superponerse y seguidamente aplica un algoritmo de programación dinámica para calcular la alineación superpuesta de puntuación máxima entre cada par restante de fragmentos. Luego los ensambla en orden de puntuación de alineación.

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

El esquema del procesado seguido se muestra en el link:

![ScreenShot](https://github.com/andreaalvarezp/Scripts_TFG/tree/main/Pipelines/pipeline_Plec.png)

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

El esquema del procesado seguido se muestra en el link:

![ScreenShot](https://github.com/andreaalvarezp/Scripts_TFG/tree/main/Pipelines/pipeline_CT.png)

A continuación se realiza una breve descripción de las utilidades que le damos a algunos de estos módulos:

#### Bowtie2 (Langmead, B., *et al,* 2012)

El módulo **Bowtie2-build** construye un índice Bowtie a partir de un conjunto de secuencias de DNA. Genera un conjunto de 6 archivos con los sufijos.1.bt2,.2.bt2, .3.bt2, .4.bt2, .rev.1.bt2 y .rev.2.bt2. que constituyen el índice. Son todo lo que se necesita para alinear las lecturas con esa referencia. Una vez creado el índice, Bowtie 2 ya no utiliza los archivos de secuencia original FASTA.

#### HMMER (Johnson, L. S. *et al,* 2010)

Herramienta que implementa métodos que utilizan modelos probabilísticos llamados modelos de perfil oculto de Markov (perfil HMM). Está diseñado para detectar homólogos remotos con la mayor sensibilidad posible, basándose en la solidez de sus modelos de probabilidad subyacentes, y es esencialmente tan rápido como BLAST.

## 2. CANDIDATOS DE SECUENCIAS DE *PLECTOSPHAERELLA* COMO MIMÉTICAS DE DOS FAMILIAS DE PÉPTIDOS FITOREGULADORES DE *ARABIDOPSIS*

### 2.1. *ssp.sh*

Este módulo se utilizó para la búsqueda de motivos SCOOP y SSP en el proteoma de *Plectosphaerella*. Uso:

```bash
$ sbatch ssp.sh $1
```

siendo ``$1`` la secuencia que se quiere utilizar como *query* para buscar en los proteomas del hongo. Posibles argumentos:
- scoop10
- scoop12
- fusarium
- verticilium
- magnaporthe
- proscoop
- SSP1 hasta SSP14

En cada uno de los archivos con el nombre del argumento deben estar alojadas la/s secuencia/s correspondientes a cada categoría. El script aplica un BLAST-P para su búsqueda y luego obtiene a partir de los resultados el archivo FASTA correspondiente que será utilizado como entrada en MEME. El módulos utilizado es:

```bash
module load BLAST+/2.9.0-gompi-2019a
```

El esquema del procesado seguido se muestra en el link:

![ScreenShot](https://github.com/andreaalvarezp/Scripts_TFG/tree/main/Pipelines/pipeline_ssp.png)

A continuación se realiza una breve descripción de las utilidades que le damos a algunos de estos módulos:

#### BLAST+ (Basic Local Alignment Search Tool, Camacho, C., *et al,* 2009)

Se trata de una herramienta que encuentra regiones de similitud local entre secuencias comparando una consulta de proteínas con una base de datos de proteínas.

#### MEME - Multiple Em for Motif Elicitation versión 5.3.3 (Bailey, T. L, *et al,* 2009)

Servidor web que proporciona un portal unificado para el descubrimiento y análisis en línea de motivos de secuencia que representan características tales como sitios de unión de ADN y dominios de interacción de proteínas. 

#### ESPript3 (Robert, X., *et al,* 2014)

Servidor web para extraer y presentar un análisis completo de la información de la estructura de la proteína primaria a cuaternaria de forma automatizada.

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

El paquete ``ComplexHeatMap`` proporciona numerosas funcionalidades para personalizar mapas de calor e incluir gráficos de anotaciones definidas por el usuario. La función ``scale()`` implementada en las matrices de datos de FPKM y expresión relativa permite que, si existen datos extremos, el efecto de aquellos datos no tan altos pero igualmente significativos no se verá diluido en la representación. 

## 3. BÚSQUEDA Y DETERMINACIÓN DE FAMILIAS DE CAZymas EN *PLECTOSPHAERELLA*

### 3.1. *CAZymas.sh*

Este módulo ejecuta múltiples tareas:
1. Ejecuta dbCAN2 para la determinación de genes que codifican para CAZymas en los proteomas y genomas de *Plectosphaerella*.
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

El esquema del procesado seguido se muestra en el link:

![ScreenShot](https://github.com/andreaalvarezp/Scripts_TFG/tree/main/Pipelines/pipeline_CAZymas.png)

A continuación se realiza una breve descripción de las utilidades que le damos a algunos de estos módulos:

#### CUPP (Barret, K., *et al,* 2020)

Herramienta de anotación funcional y agrupación no basada en alineamientos que utiliza patrones de péptidos únicos conservados para realizar agrupaciones automatizadas de proteínas y formar grupos. Es lo que se denomina una *k-mer-based tool*, donde los k-mer distintivos representan los motivos característicos de cada subfamilia y, por lo tanto, se usan para anotar nuevos genomas para CAZymas.

#### dbCAN2 (Zhang, H., *et al,* 2018)

Servidor web creado en 2012 para proporcionar un servicio público para la anotación CAZymas automatizada para genomas recién secuenciados para lo cual utiliza tres herramientas:

- **HMMER:** herramienta que implementa métodos que utilizan modelos probabilísticos llamados modelos de perfil oculto de Markov (perfil HMM).
- **DIAMOND BLAST**: herramienta de alineamiento.
- **Hotpep**: anota CAZymas mediante la búsqueda contra la biblioteca PPR (*Peptide Pattern Recognition*) de motivos peptídicos cortos conservados presente en diferentes familias CAZymas. En la biblioteca de PPR, cada familia de CAZymas tiene un conjunto de péptidos de 6-mer conservados, y Hotpep escanea nuevas proteínas para la presencia de estos péptidos con el fin de asignar las proteínas de consulta en familias CAZymas existentes.

#### SECRETOOL (Cortázar, A. R., *et al,* 2013)

Servidor web que comprende un grupo de módulos que permiten hacer predicciones de secretomas a partir de archivos de secuencias de aminoácidos. Para ello, en primer lugar, realiza un procesamiento de los datos con TargetP, SignalP y PredGPI y mezcla en un solo archivo las proteínas predichas por estos tres métodos. Seguidamente realiza una evaluación con TMHMM, método basado en cadenas de Markov. Los candidatos del paso anterior se guardan como entrada para WoLFSORT, donde las secuencias etiquetadas como “extracelulares” se retienen. Como salida se tiene una lista de IDs con las proteínas secretadas y un archivo con las secuencias relativas a esos IDs en formato FASTA. Opcionalmente, también permite predecir ortólogos y determinar dominios. 

#### IQ-TREE (Minh, B., *et al,* 2020)

Paquete de software de código abierto y ampliamente utilizado para la inferencia filogenética que utiliza el criterio de máxima verosimilitud (ML). 
