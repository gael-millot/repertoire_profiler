[//]: # "#to make links in gitlab: example with racon https://github.com/isovic/racon"
[//]: # "tricks in markdown: https://openclassrooms.com/fr/courses/1304236-redigez-en-markdown"

| usage | dependencies |
| --- | --- |
| [![Nextflow](https://img.shields.io/badge/code-Nextflow-blue?style=plastic)](https://www.nextflow.io/) | [![Dependencies: Nextflow Version](https://img.shields.io/badge/Nextflow-v22.10.3-blue?style=plastic)](https://github.com/nextflow-io/nextflow) |
| [![License: GPL-3.0](https://img.shields.io/badge/licence-GPL%20(%3E%3D3)-green?style=plastic)](https://www.gnu.org/licenses) | [![Dependencies: Singularity Version](https://img.shields.io/badge/Singularity-v3.8.0-blue?style=plastic)](https://github.com/apptainer/apptainer) |

<br /><br />
## TABLE OF CONTENTS


   - [AIM](#aim)
   - [WARNING](#warning)
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
## WARNING

Right now, only dedicated to the analysis of VDJ repertoires (corresponding to the germlines/imgt/*<SPECIES>*/vdj folder of the [IMGT database](https://www.imgt.org/IMGTrepertoire/Proteins/index.php#C)


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
DRIVE="Z"
sudo mkdir /mnt/z
sudo mount -t drvfs $DRIVE: /mnt/z
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
| **repertoires** | Folder containing the repertoires (i.e., contingency tables of the VDJ allele usage of the productive annotation by Igblast (see below). |
| **png** | Folder containing the graphs in png format. |
| **svg** | Folder containing the graphs in svg vectorial format. |
| **RData** | Folder containing, for each clonal group, objects that can be used in R to further analyze of plot the data:<br /><ul><li>db: tibble data frame resulting from the import by the alakazam::readChangeoDb() function<br /></li><li>clones: db in the airClone format<br /></li><li>trees: output of the dowser::getTrees() function using the clones object as input (igphylm tree)</li><br /><br />Also contains the all_trees.RData file that combine the trees R objects of the different files in a single trees object. |
| **seq_distance.pdf** | Distribution of the distances between the two nearest sequences (see the *nearest_distance* column in the *productive_seq.tsv* file). |
| **donuts.pdf** | donut plots showing the frequency of sequences per clonal groups, among all the productive sequences ("all", and "annotated" if metadata are provided) and among sequences used for trees (tree). |
| **repertoire.pdf** | heatmap showing the frequency of alleles used among all the productive sequences (all, non-zero and annotated sequences if metadata are provided). Non-zero means that unused alleles are removed from the heatmap (empty row or column).|
| **trees.pdf** | Phylogenic trees of the sequences that belong to a clonal group (one page per clonal group). |
| **donut_stat.tsv** | stats associated to the *donuts.pdf* file. |
| **productive_seq.tsv** | Sequences annotated by igblast (see the unproductive_seq.tsv file for sequences that failed to be productive annotated by igblast). Productive means: (1) coding region has an open reading frame, (2) no defect in the start codon, splicing sites or regulatory elements, (3) no internal stop codons, (4) an in-frame junction region. |
| **all_passed_seq.tsv** | Sequences from the *productive_seq.tsv* file, with germline clustering (clone ID), mutation load, distance and sequence nickname (annotation from the metadata file) added. Warning: the number of sequences can be lower than in the *productive_seq.tsv* file due to sequences that failed to be clone assigned (see the *non_clone_assigned_sequence.tsv* file).<br />Column description: <br /><ul><li>sequence_id: Unique query sequence identifier for the Rearrangement. Most often this will be the input sequence header or a substring thereof, but may also be a custom identifier defined by the tool in cases where query sequences have been combined in some fashion prior to alignment. When downloaded from an AIRR Data Commons repository, this will usually be a universally unique record locator for linking with other objects in the AIRR Data Model.<br /></li><li>sequence: The query nucleotide sequence. Usually, this is the unmodified input sequence, which may be reverse complemented if necessary. In some cases, this field may contain consensus sequences or other types of collapsed input sequences if these steps are performed prior to alignment.<br /></li><li>rev_comp: True if the alignment is on the opposite strand (reverse complemented) with respect to the query sequence. If True then all output data, such as alignment coordinates and sequences, are based on the reverse complement of 'sequence'.<br /></li><li>productive: True if the V(D)J sequence is predicted to be productive.<br /></li><li>v_call: V gene with allele. If referring to a known reference sequence in a database the relevant gene/allele nomenclature should be followed (e.g., IGHV4-59\*01 if using IMGT/GENE-DB).<br /></li><li>d_call: First or only D gene with allele. If referring to a known reference sequence in a database the relevant gene/allele nomenclature should be followed (e.g., IGHD3-10\*01 if using IMGT/GENE-DB).<br /></li><li>j_call: J gene with allele. If referring to a known reference sequence in a database the relevant gene/allele nomenclature should be followed (e.g., IGHJ4\*02 if using IMGT/GENE-DB).<br /></li><li>sequence_alignment: Aligned portion of query sequence, including any indel corrections or numbering spacers, such as IMGT-gaps. Typically, this will include only the V(D)J region, but that is not a requirement.<br /></li><li>germline_alignment: Assembled, aligned, full-length inferred germline sequence spanning the same region as the sequence_alignment field (typically the V(D)J region) and including the same set of corrections and spacers (if any).<br /></li><li>junction: Junction region nucleotide sequence, where the junction is defined as the CDR3 plus the two flanking conserved codons.<br /></li><li>junction_aa: Amino acid translation of the junction.<br /></li><li>v_cigar: CIGAR string for the V gene alignment. See protocol 50<br /></li><li>d_cigar: CIGAR string for the first or only D gene alignment. See protocol 50<br /></li><li>j_cigar: CIGAR string for the J gene alignment. See protocol 50<br /></li><li>stop_codon: True if the aligned sequence contains a stop codon.<br /></li><li>vj_in_frame: True if the V and J gene alignments are in-frame.<br /></li><li>locus: Gene locus (chain type). Note that this field uses a controlled vocabulary that is meant to provide a generic classification of the locus, not necessarily the correct designation according to a specific nomenclature.<br /></li><li>junction_length: Number of nucleotides in the junction sequence.<br /></li><li>np1_length: Number of nucleotides between the V gene and first D gene alignments or between the V gene and J gene alignments.<br /></li><li>np2_length: Number of nucleotides between either the first D gene and J gene alignments or the first D gene and second D gene alignments.<br /></li><li>v_sequence_start: Start position of the V gene in the query sequence (1-based closed interval).<br /></li><li>v_sequence_end: End position of the V gene in the query sequence (1-based closed interval).<br /></li><li>v_germline_start: Alignment start position in the V gene reference sequence (1-based closed interval).<br /></li><li>v_germline_end: Alignment end position in the V gene reference sequence (1-based closed interval).<br /></li><li>d_sequence_start: Start position of the first or only D gene in the query sequence. (1-based closed interval).<br /></li><li>d_sequence_end: End position of the first or only D gene in the query sequence. (1-based closed interval).<br /></li><li>d_germline_start: Alignment start position in the D gene reference sequence for the first or only D gene (1-based closed interval).<br /></li><li>d_germline_end: Alignment end position in the D gene reference sequence for the first or only D gene (1-based closed interval).<br /></li><li>j_sequence_start: Start position of the J gene in the query sequence (1-based closed interval).<br /></li><li>j_sequence_end: End position of the J gene in the query sequence (1-based closed interval).<br /></li><li>j_germline_start: Alignment start position in the J gene reference sequence (1-based closed interval).<br /></li><li>j_germline_end: Alignment end position in the J gene reference sequence (1-based closed interval).<br /></li><li>v_score: Alignment score for the V gene. See raw score<br /></li><li>v_identity: Fractional identity for the V gene alignment (proportion)<br /></li><li>v_support: V gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the V gene assignment as defined by the alignment tool.<br /></li><li>d_score: Alignment score for the first or only D gene alignment.<br /></li><li>d_identity: Fractional identity for the first or only D gene alignment.<br /></li><li>d_support: D gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the first or only D gene as defined by the alignment tool.<br /></li><li>j_score: Alignment score for the J gene alignment.<br /></li><li>j_identity: Fractional identity for the J gene alignment.<br /></li><li>j_support: J gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the J gene assignment as defined by the alignment tool.<br /></li><li>fwr1: Nucleotide sequence of the aligned FWR1 region.<br /></li><li>fwr2: Nucleotide sequence of the aligned FWR2 region.<br /></li><li>fwr3: Nucleotide sequence of the aligned FWR3 region.<br /></li><li>fwr4: Nucleotide sequence of the aligned FWR4 region.<br /></li><li>cdr1: Nucleotide sequence of the aligned CDR1 region.<br /></li><li>cdr2: Nucleotide sequence of the aligned CDR2 region.<br /></li><li>cdr3: Nucleotide sequence of the aligned CDR3 region.<br /></li><li>sequence_alignment_aa: Translation in aa of the sequence_alignment column<br /></li><li>clone_id: Clone number. A same clone_id gathers all the sequences that putatively come from a same germline cell.<br /></li><li>germline_alignment_d_mask: Because the D-segment call for B cell receptor alignments is often low confidence, the default germline format (-g dmask) places Ns in the N/P and D-segments of the junction region rather than using the D-segment assigned during reference alignment; this can be modified to generate a complete germline (-g full) or a V-segment only germline (-g vonly)<br /></li><li>germline_v_call: V germline cassette<br /></li><li>germline_d_call: D germline cassette (usually NA)<br /></li><li>germline_j_call: J germline cassette<br /></li><li>mu_count_cdr_r: number of replacement mutations in CDR1 and CDR2 of the V-segment.<br /></li><li>mu_count_cdr_s: number of silent mutations in CDR1 and CDR2 of the V-segment.<br /></li><li>mu_count_fwr_r: number of replacement mutations in FWR1, FWR2 and FWR3 of the V-segment.<br /></li><li>mu_count_fwr_s: number of silent mutations in FWR1, FWR2 and FWR3 of the V-segment.<br /></li><li>mu_count: number of replacement and silent mutations in the specified region (If frequency=TRUE and combine=TRUE)<br /></li><li>mu_freq_cdr_r: frequency of replacement mutations in CDR1 and CDR2 of the V-segment (if frequency=TRUE, R and S mutation frequencies are calculated over the number of non-N positions in the specified regions).<br /></li><li>mu_freq_cdr_s: frequency of silent mutations in CDR1 and CDR2 of the V-segment (idem).<br /></li><li>mu_freq_fwr_r: frequency of replacement mutations in FWR1, FWR2 and FWR3 of the V-segment (idem).<br /></li><li>mu_freq_fwr_s: frequency of silent mutations in FWR1, FWR2 and FWR3 of the V-segment (idem).<br /></li><li>mu_freq: frequency of replacement and silent mutations in the specified region (idem, if frequency=TRUE and combine=TRUE)<br /></li><li>dist_nearest: minimal distance from the nearest sequence using the model from the clone_model parameter (Haming by default). NA if no other sequences have same V, J and junction length or if another sequence is strictly identical (should be 0 but NA is returned)</li> |
| **unproductive_seq.tsv** | Sequences that failed productive annotations by igblast (empty file if all the sequences are productively annotated). |
| **non_clone_assigned_sequence.tsv** | Productive sequences that failed to be assigned to a clone ID by the DefineClones.py function (empty file if all the sequences are assigned). |
| **tree_clone_id.tsv** | Clonal group IDs used in the tree analysis (clonal group with at least n sequences, n being set by the nb_seq_per_clone parameter in the ig_clustering.config file). |
| **tree_dismissed_clone_id.tsv** | Clonal group IDs not used in the tree analysis (clonal group with less than n sequences, n being set by the nb_seq_per_clone parameter in the ig_clustering.config file). |
| **tree_seq.tsv** | Sequences of the *all_passed_seq.tsv* file used in the tree analysis (clonal group with at least n sequences, n being set by the nb_seq_per_clone parameter in the ig_clustering.config file). |
| **tree_dismissed_seq.tsv** | Sequences of the *all_passed_seq.tsv* file not used in the tree analysis (clonal group with less than n sequences, n being set by the nb_seq_per_clone parameter in the ig_clustering.config file). |
| **tree_seq_not_displayed.tsv** | Sequences file used in the tree analysis but not displayed in the graph, (1) because strictly identical to another sequence already in the tree and (2) because the tree_duplicate_seq parameter of the ig_clustering.config file has been set to "FALSE". |


<br /><br />
Optional elements only returned if the igblast_aa parameter is 'false' and if the input fasta are nucleotide sequences:
<br /><br />
| Ig_clustering_xxxxx folder | Description |
| --- | --- |
| **aa** | Folder containing the translation of the alignment_sequence column of the *productive_seq.tsv* file in fasta files. |
| **aligned_seq** | Folder containing the alignment_sequence column of the *productive_seq.tsv* file in fasta files. |
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


#### v8.8

quote removed from all output files


#### v8.7

xlsx2fasta.R file modified so that it now split data according to each class of the categ parameter


#### v8.6

Minor aesthetic modifications in trees and donut


#### v8.5

Bugs fixed in distToNearest


#### v8.4

- Bugs fixed in clone_assignment and get_tree with a new output file created non_clone_assigned_sequence.tsv
- Bugs fixed in tree_vizu when no metadata


#### v8.3

- Bug fixed in tree_vizu: now nb of removed seq are displayed again
- Now annotated seq are colored


#### v8.2

Bug fixed in tree_vizu: now meta are displayed again


#### v8.1

Bug fixed in get_tree


#### v8.0

- xlsx2fasta.R file modified so that fasta files have correct names
- now trees.pdf return an empty graph (but with infos) for clonal groups without trees


#### v7.1

Clean version of the xlsx2fasta.R script: can now be run on any excel file


#### v7.0

- Bug of the tree_seq_not_displayed.tsv file fixed
- Number of seq removed added in tree leafs


#### v6.11

Code debbuged. tree_seq_not_displayed.tsv remains to be debugged


#### v6.10

Repertoires improved


#### v6.9

Repertoires added


#### v6.8

- igblast_aa parameter not operational yet. igblast_aa = "true" does not work for the moment because no j data in the imgt database and no junction data are returned which block the clone_assignment process
- tree_vizu process not sensitive to cache

#### v6.7

Problem of cache fixed for distance_hist process and bug fixed for empty pdf plot in this process. It was the name seq_distance, creating a replacement of the seq_distance.pdf file


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




