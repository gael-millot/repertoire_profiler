| Usage | Requirement |
| :--- | :--- |
| [![Nextflow](https://img.shields.io/badge/code-Nextflow-blue?style=plastic)](https://www.nextflow.io/) | [![Dependencies: Nextflow Version](https://img.shields.io/badge/Nextflow-v24.10.4-blue?style=plastic)](https://github.com/nextflow-io/nextflow) |
| [![License: GPL-3.0](https://img.shields.io/badge/licence-GPL%20(%3E%3D3)-green?style=plastic)](https://www.gnu.org/licenses) | [![Dependencies: Apptainer Version](https://img.shields.io/badge/Apptainer-v1.3.5-blue?style=plastic)](https://github.com/apptainer/apptainer) |
| | [![Dependencies: Graphviz Version](https://img.shields.io/badge/Graphviz-v2.42.3-blue?style=plastic)](https://www.graphviz.org/download/) |

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
## WARNINGS

- Right now, only dedicated to the analysis of VDJ repertoires (corresponding to the germlines/imgt/*\<SPECIES\>*/vdj folder of the [IMGT database](https://www.imgt.org/IMGTrepertoire/Proteins/index.php#C).
- To make the repertoires contingency tables and heatmaps, as well as donut plots, the script currently takes the first annotation of the imgt annotation if several are presents in the *v_call*, *j_call* or *c_call* column of the *productive_seq.tsv* file.
- Only the first 100 characters of of fasta files headers are used as sequence name. And any other characters than alphanumeric are replaced by underscores.

<br /><br />
## CONTENT

| Files and folder | Description |
| :--- | :--- |
| **main.nf** | File that can be executed using a linux terminal, a MacOS terminal or Windows 10 WSL2. |
| **nextflow.config** | Parameter settings for the *main.nf* file. Users have to open this file, set the desired settings and save these modifications before execution. Of note, this configuration file is systematically saved in the reports folder (see [below](#output)) during each execution, to save the parameter settings. |
| **bin folder** | Contains files required by the *main.nf* file. |
| **Licence.txt** | Licence of the release. |


<br /><br />
## INPUT

| Required files |
| :--- |
| A folder (zipped or not) containing nucleotide fasta files, each containing a single sequence. Use table2fasta.R ([https://github.com/gael-millot/table2fasta](https://github.com/gael-millot/table2fasta)) if sequences are in a .table file. . In each fasta file, sequences can be in a several lines (\n and or \r separated). In addition, spaces and tabs can be present in the header (they will be replaced by an underscore). |
| A metadata file (optional) for adding informations in the results. |

<br />

The dataset used in the *nextflow.config* file, as an example, is available at https://zenodo.org/records/14509916/files/ig_clustering_test_1_VH.zip.

<br />

The metadata file used in the *nextflow.config* file, as an example, is available at https://zenodo.org/records/14500245/files/metadata.tsv.

<br />

Use this code to split a multi sequence fasta file into fasta files made of a single sequence:

```
FASTA_FILE="./test.fasta" # add path and name of the fasta file here
awk -v slice_size=1 -v prefix="cut" '$1 ~ /^>/{nbSeq++; currSlice=int((nbSeq-1)/slice_size)+1; myOutFile=prefix"_"currSlice".fasta"}{print $0 > myOutFile}' ${FASTA_FILE}
```

<br /><br />
## HOW TO RUN

### 1. Prerequisite

Installation of:<br />
[nextflow DSL2](https://gael-millot.github.io/protocols/docs/Protocol%20152-rev0%20DSL2.html#_Toc159933761). Please, use the version indicated above.<br />
[Graphviz](https://www.graphviz.org/download/), `sudo apt install graphviz` for Linux ubuntu.<br />
[Apptainer](https://gael-millot.github.io/protocols/docs/Protocol%20135-rev0%20APPTAINER.html#_Toc160091693).<br />
<br />

Optional installation (to avoid reccurent message) of:<br />
[Gocryptfs](https://github.com/rfjakob/gocryptfs), `sudo apt install gocryptfs` for Linux ubuntu.<br /> 

Itol key:<br />
You must freely register at https://itol.embl.de/itol_account.cgi to get your own itol key. Once registered, go to https://itol.embl.de/userInfo.cgi and click on the Toggle API access button. Then, add the key in the `phylo_tree_itolkey` parameter of the *nextflow.config* file.

<br />

### 2. Local running (personal computer)

#### 2.1. *main.nf* file in the personal computer

- Mount a server if required:

<pre>
DRIVE="Z" # change the letter to fit the correct drive
sudo mkdir /mnt/share
sudo mount -t drvfs $DRIVE: /mnt/share
</pre>

Warning: if no mounting, it is possible that nextflow does nothing, or displays a message like:
<pre>
Launching `main.nf` [loving_morse] - revision: d5aabe528b
/mnt/share/Users
</pre>

- Run the following command from where the *main.nf* and *nextflow.config* files are (example: \\wsl$\Ubuntu-20.04\home\gael):

<pre>
nextflow run main.nf -c nextflow.config
</pre>

with -c to specify the name of the config file used.

<br /><br />
#### 2.2. *main.nf* file in a public github / gitlab repository

Run the following command from where you want the results:

<pre>
nextflow run -hub pasteur gmillot/repertoire_profiler -r v1.0.0
</pre>

<br /><br />

### 3. Distant running (example with the Pasteur cluster)

#### 3.1. Pre-execution

Copy-paste this after having modified the EXEC_PATH variable:

<pre>
EXEC_PATH=$(pwd) # where the bin folder of 19583_loot is located (by default, the same path as for the main.nf file)
export CONF_BEFORE=/opt/gensoft/exe # on maestro

export JAVA_CONF=java/13.0.2
export JAVA_CONF_AFTER=bin/java # on maestro
export APP_CONF=apptainer/1.3.5
export APP_CONF_AFTER=bin/apptainer # on maestro
export GIT_CONF=git/2.39.1
export GIT_CONF_AFTER=bin/git # on maestro
export GRAPHVIZ_CONF=graphviz/2.42.3
export GRAPHVIZ_CONF_AFTER=bin/graphviz # on maestro
export GRAALVM_CONF=graalvm/ce-java17-22.3.1 # required for nextflow
export GRAALVM_CONF_AFTER=bin/graalvm # on maestro
export NEXTFLOW_CONF=nextflow/24.10.3
export NEXTFLOW_CONF_AFTER=bin/nextflow # on maestro

MODULES="${CONF_BEFORE}/${JAVA_CONF}/${JAVA_CONF_AFTER},${CONF_BEFORE}/${APP_CONF}/${APP_CONF_AFTER},${CONF_BEFORE}/${GIT_CONF}/${GIT_CONF_AFTER},${CONF_BEFORE}/${GRAPHVIZ_CONF}/${GRAPHVIZ_CONF_AFTER},${CONF_BEFORE}/${GRAALVM_CONF}/${GRAALVM_CONF_AFTER},${CONF_BEFORE}/${NEXTFLOW_CONF}/${NEXTFLOW_CONF_AFTER}"
cd ${EXEC_PATH}
chmod 755 ${EXEC_PATH}/bin/*.* # not required if no bin folder
module load ${JAVA_CONF} ${APP_CONF} ${GIT_CONF} ${GRAPHVIZ_CONF} ${GRAALVM_CONF} ${NEXTFLOW_CONF}
</pre>

<br /><br />

#### 3.2. *main.nf* file in a cluster folder

Modify the second line of the code below, and run from where the *main.nf* and *nextflow.config* files are (which has been set thanks to the EXEC_PATH variable above):

<pre>
HOME_INI=$HOME
HOME="${HELIXHOME}/repertoire_profiler/" # $HOME changed to allow the creation of .nextflow into /$HELIXHOME/repertoire_profiler/, for instance. See NFX_HOME in the nextflow software script
nextflow run main.nf -c nextflow.config # or nextflow run main.nf -c nextflow.config --modules ${MODULES} in order to have all the used module versions recorded into the report file 
HOME=$HOME_INI
</pre>

<br /><br />

#### 3.3. *main.nf* file in the public gitlab repository

Modify the first and third lines of the code below, and run (results will be where the EXEC_PATH variable has been set above):

<pre>
VERSION="v1.0"
HOME_INI=$HOME
HOME="${HELIXHOME}/repertoire_profiler/" # $HOME changed to allow the creation of .nextflow into /$HELIXHOME/repertoire_profiler/, for instance. See NFX_HOME in the nextflow software script
nextflow run -hub pasteur gmillot/repertoire_profiler -r $VERSION -c $HOME/nextflow.config
HOME=$HOME_INI
</pre>

<br /><br />
### 4. Error messages and solutions

#### Message 1
```
Unknown error accessing project `gmillot/repertoire_profiler` -- Repository may be corrupted: /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot/repertoire_profiler
```

Purge using:
<pre>
rm -rf /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot*
</pre>

#### Message 2
```
WARN: Cannot read project manifest -- Cause: Remote resource not found: https://gitlab.pasteur.fr/api/v4/projects/gmillot%2Frepertoire_profiler
```

Contact Gael Millot (distant repository is not public).

#### Message 3

```
permission denied
```

Use chmod to change the user rights. Example linked to files in the bin folder: 
```
chmod 755 bin/*.*
```

#### Message 4

```
ERROR ~ Error executing process > 'ITOL (2)'

Caused by:
  Process `ITOL (2)` terminated with an error exit status (1)
  
INFO:    underlay of /etc/localtime required more than 50 (83) bind mounts
```

Register at Itol as explained in the [Prerequisite](#1-prerequisite) section above.

<br /><br />
## OUTPUT

By default, all the results are returned in a *result* folder where the *main.nf* executed file is located (created if does not exist). This can be changed using the *out_path_ini* parameter of the *nextflow.config* file. By default, each execution produces a new folder named *repertoire_profiler_\<ID\>*, created inside the *result* folder and containing all the outputs of the execution. The name of the folder can be changed using the *result_folder_name* parameter of the *nextflow.config* file. The new name file will be followed by an \<ID\> in all cases.
<br /><br />
An example of results obtained with the dataset is present at this address: https://zenodo.org/records/15132203/files/repertoire_profiler_1743690584.zip.
<br /><br />
Complete informations are in the Protocol 144-rev0 Ig clustering - Immcantation.docx (contact Gael Millot).
<br /><br />
Mandatory elements:
<br /><br />
| repertoire_profiler_<UNIQUE_ID> folders and files | Description |
| :--- | :--- |
| **report.html** | HTML report presenting the main results of the data processing. Its content will be gradually expanded to provide a more comprehensive analysis. |
| **reports** | Folder containing all the reports of the different processes, as well as the *nextflow.config* file used and the map of the processes execution (*nf_dag.png* file). |
| **repertoires** | Folder containing the repertoires, i.e., contingency tables of the V, J and C (constant) allele and gene usage from the *clone_assigned_seq.tsv* file (see below). Warning: the script currently takes the first annotation of the imgt annotation if several are presents in the v_call, j_call and c_call columns of the *clone_assigned_seq.tsv* file. (e.g., v_call with IGKV1-39\*01,IGKV1D-39\*01), so that contingencies are identical to those from the donut frequencies, that use germline_v_call and germline_j_call columns (allele reassignment by the CreateGermlines.py tool of immcantation) |
| **figures** | Folder containing : <ul><li>**png** : Folder containing the graphs in png format. </li><li>**svg** : Folder containing the graphs in svg vectorial format.</ul> |
| **RData** | Folder containing, for each clonal group, objects that can be used in R to further analyze of plot the data:<br /><ul><li>db: tibble data frame resulting from the import by the alakazam::readChangeoDb() function<br /></li><li>clones: db in the airClone format<br /></li><li>trees: output of the dowser::getTrees() function using the clones object as input (igphylm tree)</li><br /><br />Also contains the all_trees.RData file that combine the trees R objects of the different files in a single trees object. |
| **phylo** | Folder containing alignment visualization and trees : <ul><li>**aa** : Folder containing html files displaying the alignments of the amino acid sequences in the *fasta/productive_aa* folder.<br>Each fasta file may have the following related output files:<br><ul><li>**\*_align.fasta.treefile**: Phylogenetic tree in Newick format for each sequence group in the phylo file (not constructed for groups with fewer than 4 sequences).</li><li>**\*_align.fasta_itol_url.txt**: URL to the phylogenetic tree (*\*_align.fasta.treefile*) on iTOL; tree scale is the number of mutations per position.</li></ul> </li><li>**nuc** : Folder containing html files displaying the alignments of the *"sequence_alignment"* column of the nucleotidic sequences used for germ_trees.</ul> |
| **fasta** | Folder containing all fasta files emitted in results : <ul><li>**trees_aa** : Folder containing amino-acid sequences in fasta files. Each fasta file regroups the amino-acid sequences that come from the same V\|J recombination. </li><li>**trees_nuc** : Folder containing fasta alignments of sequences from the same clonal group (e.g. with the same clone_id)</li><li>**productive_aa** : Folder containing fasta files of all producive sequences in amino-acid form. Warning: if stops are present (indicated by stars for instance), then it means that the input might come from mRNA sequencing, for which sequences can start upstream of the start codon (and up to the transcription start site), and that the column of the productive_seq.tsv used to translate is not the one to use. Please, report an issue here https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/issues</li><li>**productive_nuc** : Folder containing fasta files of all producive sequences in nucleotidic form.</li></ul> **IMPORTANT NOTE ABOUT CLONAL GROUPS** :<br>If the *clone_strategy* parameter in *nextflow.config* is set to "set", a single clonal group may contain different gene assignments. In such cases, the germline sequence included in the aligned FASTA file is selected as the most frequent germline sequence within the group.<br>If multiple germline sequences are equally frequent, the one associated with the first sequence (after sorting the group by sequence ID) is used. A warning is written to the log file whenever multiple germline sequences are detected in a group. |
| **files** | Folder containing **all the files** described below : |
| **seq_distance.pdf** | Distribution of the distances between the two nearest sequences (see the *nearest_distance* column in the *clone_assigned_seq.tsv* file). |
| **donuts.pdf** | donut plots showing the frequency of sequences per clonal groups, among:<br /><ul><li>all: all the clone-assigned sequences (*clone_assigned_seq.tsv* output file).<br /></li><li>annotated: as the "all" donut but using all the clone-assigned sequences that have been annotated using the meta_name_replacement parameter of the nextflow.config file if not "NULL".<br /></li><li>trees: all the sequences used for germline trees (*germ_tree_seq.tsv* output file).</li> |
| **\*_repertoire.pdf** | heatmap of the files from the *repertoires* folder (see above), showing the frequency of alleles and genes used among all the all clone-assigned sequences ("all"), non empty cells ("non-zero") and "annotated" sequences (if metadata are provided). Non-zero means that unused alleles are removed from the heatmap (empty row or column). Warning: to build the repertoire contingencies, the script currently takes the first annotation of the imgt annotation if several are presents in the v_call, j_call or c_call columns of the *clone_assigned_seq.tsv* file (see the *clone_assigned_seq_several_annot_igmt.tsv* file below) |
| **germ_tree.pdf** | Phylogenic trees of the sequences that belong to a clonal (supposedly germline) group made of at least n sequences, n being set by the clone_nb_seq parameter in the nextflow.config file (one page per clonal group). Warning: clonal group full names are those given by dowser::formatClones, i.e., those from germinal_v_call and germinal_j_call from the *clone_assigned_seq.tsv* file. |
| **germ_no_tree.pdf** | All the clonal groups with less than n sequences, n being set by the clone_nb_seq parameter in the nextflow.config file (one page per clonal group). Clonal group information is recapitulated in each page. |
| **donut_stat.tsv** | stats associated to the *donuts.pdf* file. |
| **igblast_unaligned_seq_name.tsv** | Names of sequences that failed to be annotated by igblast (empty file if all the sequences are annotated). |
| **igblast_aligned_seq_name.tsv** | Names of sequences annotated by igblast (more precisely by MakeDb.py igblast). If empty, generate a subsequent nextflow failure. The number lines in *igblast_unaligned_seq_name.tsv* and *igblast_aligned_seq_name.tsv* is equal to the number of submitted .fasta files. |
| **productive_seq.tsv** | Sequences annotated by igblast (see the *unproductive_seq.tsv* file for sequences that failed to be productive annotated by igblast). Productive means: (1) coding region has an open reading frame, (2) no defect in the start codon, splicing sites or regulatory elements, (3) no internal stop codons, (4) an in-frame junction region. <br /> Distance and sequence nickname (annotation from the metadata file) added. <br />Column description (from [here](https://docs.airr-community.org/en/stable/datarep/rearrangements.html#fields)): <br /><ul><li>sequence_id: Unique query sequence identifier for the Rearrangement. Most often this will be the input sequence header or a substring thereof, but may also be a custom identifier defined by the tool in cases where query sequences have been combined in some fashion prior to alignment. When downloaded from an AIRR Data Commons repository, this will usually be a universally unique record locator for linking with other objects in the AIRR Data Model.<br /></li><li>sequence: Warning: this data are not the expected sequences returned by Immcantation, which are moved into the *sequence_ini* column. In fact, the sequences are not different except a potential trimming of the 5' part that removes all the nucleotides before the fwr1 region (corresponding to the leader peptide). This should guarantees a translation of these sequences without premature stops. Trimming is performed by comparing the 5' part of the *sequence_alignment* and the *sequence_ini* sequences.<br /></li><li>rev_comp: True if the alignment is on the opposite strand (reverse complemented) with respect to the query sequence. If True then all output data, such as alignment coordinates and sequences, are based on the reverse complement of 'sequence'.<br /></li><li>productive: True if the V(D)J sequence is predicted to be productive. <br /></li><li>v_call: V cassette (gene) with allele after the star (e.g., IGHV4-59\*01 if using IMGT/GENE-DB). Sometimes, Igblast cannot distinguished between several reference sequences (cassette or allele), and are all provided, comma separated (e.g., IGHJ5\*01,IGHJ5\*02).<br /></li><li>d_call: First or only D cassette (gene) with allele after the star (e.g., IGHD3-10\*01 if using IMGT/GENE-DB). Sometimes, Igblast cannot distinguished between several reference sequences (cassette or allele), and are all provided, comma separated (e.g., IGHD3-16\*01,IGHD3-16\*02).<br /></li><li>j_call: J cassette (gene) with allele after the star (e.g., IGHJ4\*01 if using IMGT/GENE-DB). Sometimes, Igblast cannot distinguished between several reference sequences (cassette or allele), and are all provided, comma separated (e.g., IGHJ4\*02,IGHJ5\*02).<br /></li><li>sequence_alignment: Aligned portion of query sequence, including any indel corrections or numbering spacers, such as IMGT-gaps (explained [here](https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html)), in order to fit the [collier de perle represention](https://www.imgt.org/3Dstructure-DB/doc/IMGTCollier-de-Perles.shtml). Typically, this will include **only the V(D)J region**, but that is not a requirement, explaining why this sequence lacks the 3' part of the sequence indicated in the *sequence* column.<br /></li><li>germline_alignment: Assembled, aligned, full-length inferred germline sequence spanning the same region as the sequence_alignment field (typically the V(D)J region) and including the same set of corrections and IMGT-gaps spacers (if any). Thus, If well understood, this sequence is built from V, D, J reference sequences from the IMGT database that match the sequence in sequence_alignment, with IMGT-gaps spacers already present in the sequences of the database kept. The aligned parts of these V, D, J references are stitched together, and the `N` nucleotides in the germline alignment sequence are positions that belong to no templated nucleotides, added during V(D)J recombination by enzymes involved in this process, and not part of the reference. These 'N's only appear in the junction region where DJ recombination and V-DJ recombination happened and are present where VDJ reference cassettes don't match each other. This column is created by the MakeDb.py Immcantation script. Warning: sequences with the same clone_ID can have different germline_aligments since the column does not come from CreateGermlines.py. <br /></li><li>junction: Junction region nucleotide sequence, where the junction is defined as the CDR3 plus the two flanking conserved codons.<br /></li><li>junction_aa: Amino acid translation of the junction.<br /></li><li>v_cigar: CIGAR string for the V gene alignment. See protocol 50<br /></li><li>d_cigar: CIGAR string for the first or only D gene alignment. See protocol 50<br /></li><li>j_cigar: CIGAR string for the J gene alignment. See protocol 50<br /></li><li>stop_codon: True if the aligned sequence contains a stop codon.<br /></li><li>vj_in_frame: True if the V and J gene alignments are in-frame. <br /></li><li>locus: Gene locus (chain type). Note that this field uses a controlled vocabulary that is meant to provide a generic classification of the locus, not necessarily the correct designation according to a specific nomenclature.<br /></li><li>c_call: Constant region gene with allele. If referring to a known reference sequence in a database the relevant gene/allele nomenclature should be followed (e.g., IGHG1\*01 if using IMGT/GENE-DB).<br /></li><li>junction_length: Number of nucleotides in the junction sequence.<br /></li><li>np1_length: Number of nucleotides between the V gene and first D gene alignments or between the V gene and J gene alignments.<br /></li><li>np2_length: Number of nucleotides between either the first D gene and J gene alignments or the first D gene and second D gene alignments.<br /></li><li>v_sequence_start: Start position of the V gene in the query sequence (1-based closed interval).<br /></li><li>v_sequence_end: End position of the V gene in the query sequence (1-based closed interval).<br /></li><li>v_germline_start: Alignment start position in the V gene reference sequence (1-based closed interval).<br /></li><li>v_germline_end: Alignment end position in the V gene reference sequence (1-based closed interval).<br /></li><li>d_sequence_start: Start position of the first or only D gene in the query sequence. (1-based closed interval).<br /></li><li>d_sequence_end: End position of the first or only D gene in the query sequence. (1-based closed interval).<br /></li><li>d_germline_start: Alignment start position in the D gene reference sequence for the first or only D gene (1-based closed interval).<br /></li><li>d_germline_end: Alignment end position in the D gene reference sequence for the first or only D gene (1-based closed interval).<br /></li><li>j_sequence_start: Start position of the J gene in the query sequence (1-based closed interval).<br /></li><li>j_sequence_end: End position of the J gene in the query sequence (1-based closed interval).<br /></li><li>j_germline_start: Alignment start position in the J gene reference sequence (1-based closed interval). <br /></li><li>j_germline_end: Alignment end position in the J gene reference sequence (1-based closed interval).<br /></li><li>v_score: Alignment score for the V gene. See raw score<br /></li><li>v_identity: Fractional identity for the V gene alignment (proportion)<br /></li><li>v_support: V gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the V gene assignment as defined by the alignment tool.<br /></li><li>d_score: Alignment score for the first or only D gene alignment.<br /></li><li>d_identity: Fractional identity for the first or only D gene alignment.<br /></li><li>d_support: D gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the first or only D gene as defined by the alignment tool.<br /></li><li>j_score: Alignment score for the J gene alignment.<br /></li><li>j_identity: Fractional identity for the J gene alignment.<br /></li><li>j_support: J gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the J gene assignment as defined by the alignment tool.<br /></li><li>fwr1: Nucleotide sequence of the aligned FWR1 region of the query sequence (i.e., FWR1 region of the *sequence_alignment* field, which is the input sequence with IMGT gaps added).<br /></li><li>fwr2: Nucleotide sequence of the aligned FWR2 region of the query sequence (i.e., FWR2 region of the *sequence_alignment* field, which is the input sequence with IMGT gaps added).<br /></li><li>fwr3: Nucleotide sequence of the aligned FWR3 region of the query sequence (i.e., FWR3 region of the *sequence_alignment* field, which is the input sequence with IMGT gaps added).<br /></li><li>fwr4: Nucleotide sequence of the aligned FWR4 region of the query sequence (i.e., FWR4 region of the *sequence_alignment* field, which is the input sequence with IMGT gaps added).<br /></li><li>cdr1: Nucleotide sequence of the aligned CDR1 region of the query sequence (i.e., CDR1 region of the *sequence_alignment* field, which is the input sequence with IMGT gaps added).<br /></li><li>cdr2: Nucleotide sequence of the aligned CDR2 region of the query sequence (i.e., CDR2 region of the *sequence_alignment* field, which is the input sequence with IMGT gaps added).<br /></li><li>cdr3: Nucleotide sequence of the aligned CDR3 region of the query sequence (i.e., CDR3 region of the *sequence_alignment* field, which is the input sequence with IMGT gaps added).<br /></li><li>sequence_ini: the query nucleotide sequence as considered by Immcantation (this means that these sequences correspond to those expected in the *sequence* column of Immcantation outputs). Usually, this is the unmodified input sequence, which may be reverse complemented if necessary. In some cases, this field may contain consensus sequences or other types of collapsed input sequences if these steps are performed prior to alignment.<br /></li><li>is_sequence_trimmed: indicates if the sequences in the *sequence* column are trimmed or not. If FALSE, then the data in both *sequence* and *sequence_ini* are identical.<br /></li><li>sequence_aa: Translation in aa of the sequence column (e.g. query sequence potentially trimmed in the 5' part).<br /></li><li>\<OPTIONAL_COLUMN\>: reporting the data from the meta_legend parameter of the *nextflow.config* file, if ever used. Example: "KD".<br /></li><li>dist_nearest: minimal distance from the nearest sequence using the model from the clone_model parameter (Haming by default). NA if no other sequences have same V, J and junction length or if another sequence is strictly identical (should be 0 but NA is returned).<br /></li><li>v_gene: extracted from the *v_call* column but removing the allele specification after the \*.<br /></li><li>j_gene: extracted from the *j_call* column but removing the allele specification after the \*.<br /></li><li>isotype_class: extracted from the *c_call* column but indicating only the isotype.<br /></li><li>c_gene: extracted from the *c_call* column but removing the allele specification after the \*.</li> |
| **clone_assigned_seq.tsv** | Sequences from the *productive_seq.tsv* file with germline clustering (clone ID), allele reannotation (germinal_v_call and germinal_j_call columns) and mutation load added. Warning: the number of sequences (i.e., rows) can be lower than in the *productive_seq.tsv* file due to sequences that failed to be clone assigned (see the *non_clone_assigned_sequence.tsv* file).<br />Additional columns description (from [here](https://docs.airr-community.org/en/stable/datarep/rearrangements.html#fields)): <br /><ul><li>clone_id: Clone number. A same clone_id gathers all the sequences that putatively come from a same germline cell. The criteria for that is 1) same V gene, 2) same J gene, 3) same CDR3 length and 4) distance between CR3 sequences, computed using the clone_model and clone_normalize parameters from the *nextflow.config file*, less than the threshold set by the clone_distance parameter from the *nextflow.config* file<br /></li><li>germline_alignment_d_mask: as germline_alignment but with D masked (i.e., replaced by N, in the middle of the CDR3). Because the D-segment call for B cell receptor alignments is often low confidence, the default germline format (-g dmask) places Ns in the N/P and D-segments of the junction region rather than using the D-segment assigned during reference alignment; this can be modified to generate a complete germline (-g full) or a V-segment only germline (-g vonly)<br /></li><li>germline_v_call: V germline cassette<br /></li><li>germline_d_call: D germline cassette (usually NA)<br /></li><li>germline_j_call: J germline cassette<br /></li><li>germline_v_seq: nucleotide sequence of the V germline cassette<br /></li><li>germline_v_seq_no_gaps: nucleotide sequence of the V germline cassette without IMGT gaps<br /></li><li>germline_d_seq: nucleotide sequence of the D germline cassette<br /></li><li>germline_d_seq_no_gaps: nucleotide sequence of the D germline cassette without IMGT gaps<br /></li><li>germline_j_seq: nucleotide sequence of the J germline cassette<br /></li><li>germline_j_seq_no_gaps: nucleotide sequence of the J germline cassette without IMGT gaps<br /></li><li>germline_d_mask_no_gaps: *germline_alignment_d_mask* column without IMGT gaps<br /></li><li>germline_d_mask_aa_no_gaps: translation of the *germline_d_mask_no_gaps* column into amino-acids<br /></li><li>mu_count_cdr_r: number of replacement mutations in CDR1 and CDR2 of the V-segment (comparing column *sequence_alignment* and *column germline_alignment_d_mask*, see https://shazam.readthedocs.io/en/stable/topics/observedMutations/#value).<br /></li><li>mu_count_cdr_s: number of silent mutations in CDR1 and CDR2 of the V-segment.<br /></li><li>mu_count_fwr_r: number of replacement mutations in FWR1, FWR2 and FWR3 of the V-segment.<br /></li><li>mu_count_fwr_s: number of silent mutations in FWR1, FWR2 and FWR3 of the V-segment.<br /></li><li>mu_count: number of replacement and silent mutations (sum of the previous columns).<br /></li><li>mu_freq_cdr_r: frequency of replacement mutations in CDR1 and CDR2 of the V-segment (if frequency=TRUE, R and S mutation frequencies are calculated over the number of non-N positions in the specified regions).<br /></li><li>mu_freq_cdr_s: frequency of silent mutations in CDR1 and CDR2 of the V-segment (idem).<br /></li><li>mu_freq_fwr_r: frequency of replacement mutations in FWR1, FWR2 and FWR3 of the V-segment (idem).<br /></li><li>mu_freq_fwr_s: frequency of silent mutations in FWR1, FWR2 and FWR3 of the V-segment (idem).<br /></li><li>mu_freq: frequency of replacement and silent mutations (sum of the previous columns).</li> |
| **clone_assigned_seq_several_annot_igmt.tsv** | Sequences from the *clone_assigned_seq.tsv* file with several annotation in the v_call and j_call columns, because several alleles had the best alignment using imgt blast (same ) |
| **unproductive_seq.tsv** | Sequences that failed productive annotations by igblast (empty file if all the sequences are productively annotated). |
| **non_clone_assigned_sequence.tsv** | Productive sequences that failed to be assigned to a clone ID by the DefineClones.py function (empty file if all the sequences are assigned). See https://changeo.readthedocs.io/en/latest/methods/clustering.html but failure reasons are not explained.|
| **germ_tree_clone_id.tsv** | Clonal group IDs used in the germline tree analysis (clonal group with at least n sequences, n being set by the clone_nb_seq parameter of the *nextflow.config* file). |
| **germ_tree_dismissed_clone_id.tsv** | Clonal group IDs not used in the germline tree analysis (clonal group with less than n sequences, n being set by the clone_nb_seq parameter of the *nextflow.config* file). |
| **germ_tree_seq.tsv** | Sequences of the *clone_assigned_seq.tsv* file used in the germline tree analysis (clonal group with at least n sequences, n being set by the clone_nb_seq parameter of the *nextflow.config* file). Additionnal germline_[*vdj*]_gene columns. |
| **germ_tree_dismissed_seq.tsv** | Sequences of the *clone_assigned_seq.tsv* file not used in the germline tree analysis (clonal group with less than n sequences, n being set by the clone_nb_seq parameter of the *nextflow.config* file). |
| **germ_tree_dup_seq_not_displayed.tsv** | Sequences file used in the germline tree analysis but not displayed in the graph, (1) because strictly identical to another sequence already in the tree and (2) because the tree_duplicate_seq parameter of the *nextflow.config* file has been set to "FALSE". |


<br /><br />
Optional elements only returned if the igblast_aa parameter is 'false' and if the input fasta are nucleotide sequences:
<br /><br />
| repertoire_profiler_<UNIQUE_ID> folders and files | Description |
| :--- | :--- |
| **aa** | Folder containing the translation of the alignment_sequence column of the *productive_seq.tsv* file in fasta files. |
| **aligned_seq** | Folder containing the alignment_sequence column of the *productive_seq.tsv* file in fasta files. |
| **aa.tsv** |  File containing all the translation of the alignment_sequence column of the *productive_seq.tsv* file. |


<br /><br />
## VERSIONS


The different releases are tagged [here](https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/tags).

<br /><br />
## LICENCE


This package of scripts can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchandability or fitness for a particular purpose.
See the GNU General Public License for more details at https://www.gnu.org/licenses or in the Licence.txt attached file.


<br /><br />
## CITATION


Version V10.3:

[Dejoux A, Zhu Q, Ganneau C, Goff OR, Godon O, Lemaitre J, Relouzat F, Huetz F, Sokal A, Vandenberghe A, Pecalvel C, Hunault L, Derenne T, Gillis CM, Iannascoli B, Wang Y, Rose T, Mertens C, Nicaise-Roland P; NASA Study Group; England P, Mahévas M, de Chaisemartin L, Le Grand R, Letscher H, Saul F, Pissis C, Haouz A, Reber LL, Chappert P, Jönsson F, Ebo DG, Millot GA, Bay S, Chollet-Martin S, Gouel-Chéron A, Bruhns P. Rocuronium-specific antibodies drive perioperative anaphylaxis but can also function as reversal agents in preclinical models. Sci Transl Med. 2024 Sep 11;16(764):eado4463. doi: 10.1126/scitranslmed.ado4463. Epub 2024 Sep 11. PMID: 39259810.](https://www.science.org/doi/10.1126/scitranslmed.ado4463)


<br /><br />
## CREDITS

[Pascal Chappert](https://www.institut-necker-enfants-malades.fr/index.php?menu=team&rubric=teamtabs&idfac=mahevas#chart), INSERM U1151 Institut Necker Enfants Malades, Paris, France

[Frédéric Lemoine](), Institut Pasteur, Université Paris Cité, Bioinformatics and Biostatistics Hub, 75015 Paris, France

[Chloé Taurel](), University of Technology of Compiègne, Compiègne, France

[Gael A. Millot](https://gitlab.pasteur.fr/gmillot), Institut Pasteur, Université Paris Cité, Bioinformatics and Biostatistics Hub, 75015 Paris, France

<br /><br />
## ACKNOWLEDGEMENTS


The developers & maintainers of the mentioned softwares and packages, including:

- [R](https://www.r-project.org/)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [immacantation](https://immcantation.readthedocs.io/en/stable/)
- [ggtree](https://yulab-smu.top/treedata-book/)
- [Nextflow](https://www.nextflow.io/)
- [Apptainer](https://apptainer.org/)
- [Docker](https://www.docker.com/)
- [Gitlab](https://about.gitlab.com/)
- [Bash](https://www.gnu.org/software/bash/)
- [Ubuntu](https://ubuntu.com/)

Special acknowledgement to [Kenneth Hoehn](https://medicine.yale.edu/profile/kenneth-hoehn/), Yale School of Medicine, New Haven, CT, USA

<br /><br />
## WHAT'S NEW IN

#### v19.2

-  Bug fixed

#### v19.1

-  Bug fixed and README completed

#### v19.0

-  Now sequences are trimmed in 5' for the leader peptide if exists, and two new columns in productive_seq.tsv sequence_ini is_sequence_trimmed"

#### v18.6

- igblast_aa parameter removed from the nextflow.config file

#### v18.5

- Many things improved or added

#### v18.4

- goalign version updated for correct running in *nextflow.config*

#### v18.3

- Trees and alignments made on the *sequence* column (query sequence) instead of the *sequence_alignment* column
- germ_tree result files removed
- Folder reorganization
- Columns added to *clone_assigned_seq.tsv* result file
- New *.gff* files for region coordinates in nucleotidic alignment visualization (in progess)
- Bugs fixed for meta_path = "NULL" in *nextflow.config*

#### v18.2

- Bug fixes for gene names containing "/" character

#### v18.1

- tsv files truncated to 10 lines in html report
- Bug fixes for HPC use

#### v18.0

- Tool adapted to any cluster job scheduler
- *phylo_tree_heavy* parameter removed from config, chain type is now deduced from igblast_ref files
- Processes *Reformat*, *Align*, *DefineGroups*, *NbSequences*, *Tree*, *ProcessMeta* and *ITOL* are no longer restricted to heavy chain analysis
- cdr3 alignments removed (in development for later release)
- tsv files better displayed in html report

#### v17.0

- New config parameter added for ITOL subscription
- New alignment folders displayed in html for cdr3 and aminoacid sequences

#### v16.5

- Bug fix in print_report for larger data input

#### v16.4

- Display improved in repertoires
- Identical sequences must be kept in trees by default (germ_tree_duplicate_seq parameter in nextflow.config)

#### v16.3

- Precisions on K and L loci incorporated in donuts and repertoires for light chain analysis

#### v16.2

- HTML report now includes donuts, repertoires, statistics on successfull sequence outputs following several processes, and clonal  prediction
- Bugs regarding the study of K and L loci on light chains fixed
- Thorough checking of *ig_ref_files* format

#### v16.0

- HTML report produced to present the main results in a more accessible way (work in progress)
- All result files located in folders (new *"file"* folder)
- *productive_seq.tsv* updated to contain translation and name_replacement information


#### v15.0

- Now repertoire files (.tsv and .pdf) deal with the constant chain.


#### v14.3

- Isotype class and subclass columns added in the output clone_assigned_seq.tsv file.


#### v14.2

- Bug fixed in the metadata display in trees.
- Tree display improved.


#### v14.1

- Bug fixed in the metadata display in trees.
- New metadata set.


#### v14.0

- New meta_seq_names parameter in the nextflow.config file, that help to set the name of the column containing the fasta sequence names.


#### v13.0

- Immcantation upgraded to 4.5.0: now the constant region allele is provided


#### v12.0

- Deal with nextflow v24.10.4, test GString in addition of String for path


#### v11.9

- Internal error messages improved


#### v11.8

- Error messages improved so that now nextflow does not display all subsequent error messages because of empty file when execution is stopped by an error


#### v11.7

- Bug solved in the translation process: now parallelization works even if empty input file


#### v11.6

- Fasta files can have tabs in the header


#### v11.5

- nextflow.config file improved for users
- README file updated
- Title of donut improved
- Zenodo output link updated


#### v11.4

- Zenodo input link updated


#### v11.3

- Now a zip folder is accepted as input
- Explanations in nextflow.config are improved


#### v11.2

- Donut and Tree plots slightly improved


#### v11.1

- Option comment_char = '' of read.table() for all .R files with read.table() in bin


#### v11.0

- Tree plots improved


#### v10.4

- Donut plot legend corrected


#### v10.3

- Spaces removed in file names
- Spaces removed in the first line of each fasta file
- .fas fasta extension allowed


#### v10.2

xlsx2fasta.R file: error fixed for the fasta by categ


#### v10.1

xlsx2fasta.R file improved to take into account the problem of NA


#### v10.0

xlsx2fasta.R file strongly improved to deal with empty sequences


#### v9.9

Bug fixed in xlsx2fasta.R


#### v9.8

Error fixed in the nextflow.config file about cute path


#### v9.7

Cute folder updated to 12.8 to fix the palette bug in fun_gg_donut()


#### v9.6

Bug in the metadata_check process fixed


#### v9.5

README improved so that now dataset and results are in zenodo


#### v9.4

bug fixed


#### v9.3

important check added for the metadata file


#### v9.2

- repertoires now ok with a mix of IGK and IGL sequences
- first annotation taken if several v or j allele annotation by imgt for the repertoires


#### v9.1

bug corrected for the meta_name_replacement parameter


#### v9.0

repertoire_profiler.nf and .config name changed so that now can be run from gitlab


#### v8.14

README file improved, that clarify the differences between sequence_alignment and germline_\* sequences


#### v8.13

README file improved, that clarify if results are from the productive seq or all passed seq


#### v8.12

Check added for the metadata file of the meta_path parameter. But check the content of the first column of this file remains to be added


#### v8.11

bugs fixed


#### v8.10

bugs fixed


#### v8.9

ig_clustering name replaced everywhere by repertoire_profiler


#### v8.8

Quotes removed from all output files


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



