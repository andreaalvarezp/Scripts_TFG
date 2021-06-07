# README: Búsqueda de motivos genómicos con posible relevancia biológica en tres aislados de hongos del género *Plectosphaerella*

**Grado en Biotecnología, Universidad Politécnica de Madrid**

**Autor: Andrea Álvarez Pérez**

**Curso 2020-2021**

Este Github recoge los scripts bioinformáticos utilizados para la realización del **Trabajo de Fin de Grado de Biotecnología** de la alumna Andrea Álvarez Pérez. Cada uno de ellos está destinado a una tarea concreta y todos los scripts, tanto en R como en bash, están comentados con la función y el objetivo de cada una de las órdenes y módulos ejecutados.

De igual manera, en este archivo se explica de manera concisa el uso de cada uno de los script y se clasifican según la sección del Trabajo de Fin de Grado a la que pertenecen.

## 1. Búsqueda y determinación de micovirus en *Plectosphaerella*

### *Plec.sh*

Este módulo aplica un procesa bioinformático a los datos de RNA-seq brutos de *Plectosphaerella* que permite la identificación de secuencias de micovirus en los mismos. Uso:

```bash
sbatch Plec.sh $1
```

siendo ``$1`` la cepa de *Plectosphaerella* a analizar:
- PcBMM
- Pc2127
- P0831





