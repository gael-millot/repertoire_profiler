[//]: # "#to make links in gitlab: example with racon https://github.com/isovic/racon"
[//]: # "tricks in markdown: https://openclassrooms.com/fr/courses/1304236-redigez-en-markdown"

| usage | dependencies |
| --- | --- |
| [![Nextflow](https://img.shields.io/badge/code-Nextflow-blue?style=plastic)](https://www.nextflow.io/) | [![Dependencies: Nextflow Version](https://img.shields.io/badge/Nextflow-v22.10.3-blue?style=plastic)](https://github.com/nextflow-io/nextflow) |
| [![License: GPL-3.0](https://img.shields.io/badge/licence-GPL%20(%3E%3D3)-green?style=plastic)](https://www.gnu.org/licenses) | [![Dependencies: Singularity Version](https://img.shields.io/badge/Singularity-v3.8.0-blue?style=plastic)](https://github.com/apptainer/apptainer) |

<br /><br />
## TABLE OF CONTENTS


   - [AIM](#aim)
   - [CONTENT](#content)
   - [INPUT](#input)
   - [HOW TO RUN](#how-to-run)
   - [OUTPUT](#output)
   - [VERSIONS](#versions)
   - [LICENCE](#licence)
   - [CITATION](#citation)
   - [CREDITS](#credits)
   - [ACKNOWLEDGEMENTS](#Acknowledgements)
   - [WHAT'S NEW IN](#what's-new-in)

<br /><br />
## AIM

- Annotation of mRNA sequencing of the Immunoglobuline Heavy or Light variable region (fasta file sequences).
- Clustering of the annotated sequences into clonal groups (same germline origin).
- Tree visualization of the clonal groups.


<br /><br />
## CONTENT

| Ig_clustering folder | Description |
| --- | --- |
| **ig_clustering.nf** | File that can be executed using a linux terminal, a MacOS terminal or Windows 10 WSL2. |
| **ig_clustering.config** | Parameter settings for the ig_clustering.nf file. Users have to open this file, set the desired settings and save these modifications before execution. |
| **xlsx2fasta.R** | Accessory file that creates all the fasta files from a .xlsx file. To use it, 1) open the file, 2) complete the "Parameters that need to be set by the user" section, 3) save the modifications and 4) run the file in R. |
| **dataset** | Folder containing some datasets (batch of fasta files) than can be used as examples. |
| **example_of_results** | Folder containing examples of result obtained with the dataset. See the OUTPUT section for the description of the folder and files. |


<br /><br />
## INPUT

A folder made of fasta files, each containing a single sequence.

Use the xlsx2fasta.R script if sequences are in a .xlsx file (see the sequences.xlsx file in the dataset folder as an example).

Use this code to split a multi sequence fasta file into fasta files made of a single sequence:

<pre>
FASTA_FILE="./test.fasta" # add path and name of the fasta file here
awk -v slice_size=1 -v prefix="cut" '$1 ~ /^>/{nbSeq++; currSlice=int((nbSeq-1)/slice_size)+1; myOutFile=prefix"_"currSlice".fasta"}{print $0 > myOutFile}' ${FASTA_FILE}
</pre>

<br /><br />
## HOW TO RUN

### 1. Prerequisite

Installation of:<br />
[nextflow DSL2](https://github.com/nextflow-io/nextflow)<br />
[Graphviz](https://www.graphviz.org/download/), `sudo apt install graphviz` for Linux ubuntu<br />
[Singularity/apptainer](https://github.com/apptainer/apptainer)<br />

<br /><br />
### 2. Local running (personal computer)


#### 2.1. ig_clustering.nf file in the personal computer

- Mount a server if required:

<pre>
DRIVE="C"
sudo mkdir /mnt/c
sudo mount -t drvfs $DRIVE: /mnt/c
</pre>

Warning: if no mounting, it is possible that nextflow does nothing, or displays a message like:
<pre>
Launching `ig_clustering.nf` [loving_morse] - revision: d5aabe528b
/mnt/share/Users
</pre>

- Run the following command from where the ig_clustering.nf and ig_clustering.config files are (example: \\wsl$\Ubuntu-20.04\home\gael):

<pre>
nextflow run ig_clustering.nf -c ig_clustering.config
</pre>

with -c to specify the name of the config file used.

<br /><br />
#### 2.3. ig_clustering.nf file in the public gitlab repository

Run the following command from where you want the results:

<pre>
nextflow run -hub pasteur gmillot/ig_clustering -r v1.0.0
</pre>

<br /><br />
### 3. Distant running (example with the Pasteur cluster)

#### 3.1. Pre-execution

Copy-paste this after having modified the EXEC_PATH variable:

<pre>
EXEC_PATH="/pasteur/zeus/projets/p01/BioIT/gmillot/ig_clustering" # where the bin folder of the ig_clustering.nf script is located
export CONF_BEFORE=/opt/gensoft/exe # on maestro

export JAVA_CONF=java/13.0.2
export JAVA_CONF_AFTER=bin/java # on maestro
export SINGU_CONF=apptainer/1.1.5
export SINGU_CONF_AFTER=bin/singularity # on maestro
export GIT_CONF=git/2.39.1
export GIT_CONF_AFTER=bin/git # on maestro
export GRAPHVIZ_CONF=graphviz/2.42.3
export GRAPHVIZ_CONF_AFTER=bin/graphviz # on maestro

MODULES="${CONF_BEFORE}/${JAVA_CONF}/${JAVA_CONF_AFTER},${CONF_BEFORE}/${SINGU_CONF}/${SINGU_CONF_AFTER},${CONF_BEFORE}/${GIT_CONF}/${GIT_CONF_AFTER}/${GRAPHVIZ_CONF}/${GRAPHVIZ_CONF_AFTER}"
cd ${EXEC_PATH}
# chmod 755 ${EXEC_PATH}/bin/*.* # not required if no bin folder
module load ${JAVA_CONF} ${SINGU_CONF} ${GIT_CONF} ${GRAPHVIZ_CONF}
</pre>

<br /><br />
#### 3.2. ig_clustering.nf file in a cluster folder

Modify the second line of the code below, and run from where the ig_clustering.nf and ig_clustering.config files are (which has been set thanks to the EXEC_PATH variable above):

<pre>
HOME_INI=$HOME
HOME="${ZEUSHOME}/ig_clustering/" # $HOME changed to allow the creation of .nextflow into /$ZEUSHOME/ig_clustering/, for instance. See NFX_HOME in the nextflow software script
trap '' SIGINT
nextflow run --modules ${MODULES} ig_clustering.nf -c ig_clustering.config
HOME=$HOME_INI
trap SIGINT
</pre>

<br /><br />
#### 3.3. ig_clustering.nf file in the public gitlab repository

Modify the first and third lines of the code below, and run (results will be where the EXEC_PATH variable has been set above):

<pre>
VERSION="v1.0"
HOME_INI=$HOME
HOME="${ZEUSHOME}/ig_clustering/" # $HOME changed to allow the creation of .nextflow into /$ZEUSHOME/ig_clustering/, for instance. See NFX_HOME in the nextflow software script
trap '' SIGINT
nextflow run --modules ${MODULES} -hub pasteur gmillot/ig_clustering -r $VERSION -c $HOME/ig_clustering.config
HOME=$HOME_INI
trap SIGINT
</pre>

<br /><br />
### 4. Error messages and solutions

#### Message 1
```
Unknown error accessing project `gmillot/ig_clustering` -- Repository may be corrupted: /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot/ig_clustering
```

Purge using:
<pre>
rm -rf /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot*
</pre>

#### Message 2
```
WARN: Cannot read project manifest -- Cause: Remote resource not found: https://gitlab.pasteur.fr/api/v4/projects/gmillot%2Fig_clustering
```

Contact Gael Millot (distant repository is not public).

#### Message 3

```
permission denied
```

Use chmod to change the user rights.


<br /><br />
## OUTPUT


Results are present in a *Ig_clustering_xxxxx* folder, inside the nextflow *result* folder.
<br /><br />
Complete informations are in the Protocol 144-rev0 Ig clustering - Immcantation.docx (contact Gael Millot).
<br /><br />


Mandatory elements:
<br /><br />
| Ig_clustering_xxxxx folder | Description |
| --- | --- |
| **reports** | Folder containing all the reports of the different processes, including the *ig_clustering.config* file used. |
| **png** | Folder containing the graphs in png format. |
| **svg** | Folder containing the graphs in svg vectorial format. |
| **RData** | Folder containing, for each clonal group, objects that can be used in R to further analyze of plot the data:<br /><ul><li>db: tibble data frame resulting from the import by the alakazam::readChangeoDb() function<br /></li><li>clones: db in the airClone format<br /></li><li>trees: output of the dowser::getTrees() function using the clones object as input (igphylm tree)</li><br /><br />Also contains the all_trees.RData file that combine the trees R objects of the different files in a single trees object. |
| **nearest_distance.tsv** | File showing the haming distances, in the last "dist_nearest" column, used in the hamming_distance.pdf plot |
| **hamming_distance.pdf** | Distribution of the Haming distances between the two nearest sequences (see the *nearest_distance.tsv* file). |
| **productive_seq.tsv** | Sequences annotated by igblast, with germline clustering and mutation load added (see the unproductive_seq.tsv file for sequences that failed to be annotated by igblast). |
| **unproductive_seq.tsv** | Sequences that failed annotations by igblast (empty file if all the sequences are annotated). |
| **tree_clone_id.tsv** | Clonal group IDs used in the tree analysis (clonal group with at least n sequences, n being set by the nb_seq_per_clone parameter in the ig_clustering.config file). |
| **tree_dismissed_clone_id.tsv** | Clonal group IDs not used in the tree analysis (clonal group with less than n sequences, n being set by the nb_seq_per_clone parameter in the ig_clustering.config file). |
| **tree_seq.tsv** | Sequences of the *productive_seq.tsv* file used in the tree analysis (clonal group with at least n sequences, n being set by the nb_seq_per_clone parameter in the ig_clustering.config file). |
| **tree_dismissed_seq.tsv** | Sequences of the *productive_seq.tsv* file not used in the tree analysis (clonal group with less than n sequences, n being set by the nb_seq_per_clone parameter in the ig_clustering.config file). |
| **tree_seq_not_displayed.tsv** | Sequences file used in the tree analysis but not displayed in the graph, (1) because strictly identical to another sequence already in the tree and (2) because the tree_duplicate_seq parameter of the ig_clustering.config file has been set to "FALSE". |
| **trees.pdf** | Phylogenic trees of the sequences that belong to a clonal group (one page per clonal group). |
| **donuts.pdf** | Frequency of sequences per clonal groups, among sequences used for trees (tree), among all the productive sequences (all) and among functional sequences (functional). |

Optional elements only returned if the igblast_aa parameter is 'false' and if the input fasta are nucleotide sequences:
<br /><br />
| Ig_clustering_xxxxx folder | Description |
| --- | --- |
| **aligned_seq** | Folder containing the alignment_sequence column of the *productive_seq.tsv* file in fasta files. |
| **aa** | Folder containing the translation of the alignment_sequence column of the *productive_seq.tsv* file in fasta files. |
| **aa.tsv** |  File containing all the translation of the alignment_sequence column of the *productive_seq.tsv* file. |


<br /><br />
## VERSIONS


The different releases are tagged [here](https://gitlab.pasteur.fr/gmillot/ig_clustering/-/tags)

<br /><br />
## LICENCE


This package of scripts can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchandability or fitness for a particular purpose.
See the GNU General Public License for more details at https://www.gnu.org/licenses.

<br /><br />
## CITATION


Not yet published


<br /><br />
## CREDITS

[Pascal Chappert](https://www.institut-necker-enfants-malades.fr/index.php?menu=team&rubric=teamtabs&idfac=mahevas#chart), INSERM U1151 Institut Necker Enfants Malades, Paris, France

[Gael A. Millot](https://gitlab.pasteur.fr/gmillot), Bioinformatics and Biostatistics Hub, Institut Pasteur, Paris, France

<br /><br />
## ACKNOWLEDGEMENTS


The developers & maintainers of the mentioned softwares and packages, including:

- [The R immcantation solution](https://immcantation.readthedocs.io/en/stable/)
- [Bash](https://www.gnu.org/software/bash/)
- [R](https://www.r-project.org/)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [ggtree](https://yulab-smu.top/treedata-book/)
- [Nextflow](https://www.nextflow.io/)
- [Singularity/Apptainer](https://apptainer.org/)
- [Docker](https://www.docker.com/)
- [Gitlab](https://about.gitlab.com/)

Special acknowledgement to [Kenneth Hoehn](https://medicine.yale.edu/profile/kenneth-hoehn/), Yale School of Medicine, New Haven, CT, USA

<br /><br />
## WHAT'S NEW IN


#### v6.6

Names of metadata now systematically present in trees, even when identical sequences are removed, and bug fixed for empty pdf plot


#### v6.5

Code secured and bug fixed for NULL metadata file path


#### v6.4

Distances added in the returned productive_seq.tsv, now 3 histograms to help to set the distance threshold, many bugs fixed


#### v6.3

AA sequences added. Translation of the alignment_sequence column is added in the returned productive_seq.tsv and new aa.tsv file


#### v6.2

Supp donut added regarding clonal groups with functional annotation


#### v6.1

New parameters for the donut


#### v6.0

Distance histogram added
Empty graphs added
New parameters for the donut


#### v5.1

Metadata info added in donut plot


#### v5.0

Bug solved when the tree_meta_path_names parameter is a categorical column
New tree_meta_name_replacement parameter

#### v4.4

Bug solved in seq_not_displayed.tsv file


#### v4.3

igblast_aa parameter removed from the .config file
bug solved in tree_vizu


#### v4.2

Donut charts grouped in a single pdf


#### v4.1

seq_not_displayed.tsv file added to better understand the absence of seq in trees<br />
V and J alleles added in the donut legend


#### v4.0

tree_meta_path modified so that it can now be a 'NULL' path<br />
V and J alleles added in tree titles<br />
Duplicated sequences can be removed or not from trees


#### v3.5

Bug solved


#### v3.4

Two Donut charts added


#### v3.3

Empty channel solved<br />
Donut chart added


#### v3.2

README file updated for localization of the igblast database


#### v3.1

README file updated


#### v3.0

First version that provides trees of clonal groups


#### v2.1

xlsx2fasta.R file added | Nicer tree representation added


#### v2.0

Conversion into DSL2 ok


#### v1.0

First DSL1 version that works




