# README: Búsqueda de motivos genómicos con posible relevancia biológica en tres aislados de hongos del género *Plectosphaerella*

**Grado en Biotecnología, Universidad Politécnica de Madrid**

**Autor: Andrea Álvarez Pérez**

**Curso 2020-2021**

Este Github recoge los scripts bioinformáticos utilizados para la realización del **Trabajo de Fin de Grado de Biotecnología** de la alumna Andrea Álvarez Pérez. Cada uno de ellos está destinado a una tarea concreta y todos los scripts, tanto en R como en bash, están comentados con la función y el objetivo de cada una de las órdenes y módulos ejecutados.

De igual manera, en este archivo se explica de manera concisa el uso de cada uno de los script y se clasifican según la sección del Trabajo de Fin de Grado a la que pertenecen. 

**Importante**: si se quiere hacer uso de los scripts es obligatorio modificar los *paths* establecidos en el primer bloque de cada script y redefinir las variables de manera personalizada para cada usuario.

## 1. Búsqueda y determinación de micovirus en *Plectosphaerella*

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

El uso y los parámetros de cada módulo se especifican dentro del script.

### 1.2. ct.sh

Este módulo es idéntico al anterior, pero posee las variables redefinidas para utilizarlas con datos del hongo *Colletotrichum tofieldiae*. Sirve para la búsqueda e identificación de micovirus en datos de RNAseq. Uso:

```bash
$ sbatch ct.sh $1
```

siendo ``$1`` las condiciones de las muestras de *Colletotrichum tofieldiae*. Posibles argumentos:
- *invitro*
- plusP
- minusP

