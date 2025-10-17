| Usage | Requirement |
| :--- | :--- |
| [![Nextflow](https://img.shields.io/badge/code-Nextflow-blue?style=plastic)](https://www.nextflow.io/) | [![Dependencies: Nextflow Version](https://img.shields.io/badge/Nextflow-v24.10.4-blue?style=plastic)](https://github.com/nextflow-io/nextflow) |
| [![License: GPL-3.0](https://img.shields.io/badge/licence-GPL%20(%3E%3D3)-green?style=plastic)](https://www.gnu.org/licenses) | [![Dependencies: Apptainer Version](https://img.shields.io/badge/Apptainer-v1.3.5-blue?style=plastic)](https://github.com/apptainer/apptainer) |
| | [![Dependencies: Graphviz Version](https://img.shields.io/badge/Graphviz-v2.42.3-blue?style=plastic)](https://www.graphviz.org/download/) |

<br><br>
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

<br><br>
## AIM

- Annotation of mRNA sequencing of the Immunoglobuline Heavy or Light variable region (fasta file sequences).
- Clustering of the annotated sequences into clonal groups (same germline origin).
- Tree visualization of the clonal groups.

<br><br>
## WARNINGS

- Right now, only dedicated to the analysis of VDJ repertoires (corresponding to the germlines/imgt/*\<SPECIES\>*/vdj folder of the [IMGT database](https://www.imgt.org/IMGTrepertoire/Proteins/index.php#C).
- To make the repertoires contingency tables and heatmaps, as well as donut plots, the script currently takes the first annotation of the imgt annotation if several are presents in the *v_call*, *j_call* or *c_call* column of the *productive_seq.tsv* file.
- Only the first 100 characters of of fasta files headers are used as sequence name. And any other characters than alphanumeric are replaced by underscores.

<br><br>
## CONTENT

| Files and folder | Description |
| :--- | :--- |
| **main.nf** | File that can be executed using a linux terminal, a MacOS terminal or Windows 10 WSL2. |
| **nextflow.config** | Parameter settings for the *main.nf* file. Users have to open this file, set the desired settings and save these modifications before execution. Of note, this configuration file is systematically saved in the reports folder (see [below](#output)) during each execution, to save the parameter settings. |
| **bin folder** | Contains files required by the *main.nf* file. |
| **Licence.txt** | Licence of the release. |


<br><br>
## INPUT

| Required files |
| :--- |
| A folder (zipped or not) containing nucleotide fasta files, each containing a single sequence.<br>Use table2fasta.R ([https://github.com/gael-millot/table2fasta](https://github.com/gael-millot/table2fasta)) if sequences are in a .table file.<br>In each fasta file, the sequence can be split into several lines (\n and or \r separated). In addition, spaces and tabs can be present in the header (they will be replaced by an underscore). |
| A metadata file (optional) for adding informations in the results. |

<br>

The dataset used in the *nextflow.config* file, as an example, is available at https://zenodo.org/records/14509916/files/ig_clustering_test_1_VH.zip.

<br>

The metadata file used in the *nextflow.config* file, as an example, is available at https://zenodo.org/records/14500245/files/metadata.tsv.

<br>

Use this code to split a multi sequence fasta file into fasta files made of a single sequence:

```
FASTA_FILE="./test.fasta" # add path and name of the fasta file here
awk -v slice_size=1 -v prefix="cut" '$1 ~ /^>/{nbSeq++; currSlice=int((nbSeq-1)/slice_size)+1; myOutFile=prefix"_"currSlice".fasta"}{print $0 > myOutFile}' ${FASTA_FILE}
```

<br><br>
## HOW TO RUN

### 1. Prerequisite

Installation of:<br>
[nextflow DSL2](file:///C:/Users/gmillot/Documents/Git_projects/protocols/docs/Protocol%20152-rev0%20DSL2.htm#_Toc208504071). Please, use the version indicated above.<br>
[Graphviz](https://www.graphviz.org/download/), `sudo apt install graphviz` for Linux ubuntu.<br>
[Apptainer](https://gael-millot.github.io/protocols/docs/Protocol%20135-rev0%20APPTAINER.html#_Toc160091693).<br>
<br>

Optional installation (to avoid reccurent message) of:<br>
[Gocryptfs](https://github.com/rfjakob/gocryptfs), `sudo apt install gocryptfs` for Linux ubuntu.<br> 

Itol key:<br>
If you need sequence phylogenic trees in the output, you can freely register at https://itol.embl.de/itol_account.cgi to get your own itol key. Once registered, go to https://itol.embl.de/userInfo.cgi and click on the Toggle API access button. Then, add the key in the `phylo_tree_itolkey` parameter of the *nextflow.config* file, and set the `phylo_tree_itol_subscription` parameter to `TRUE`. If you experience problem with registration, set the `phylo_tree_itol_subscription` parameter to `FALSE`. The html output file explains how to see the trees without ITOL key.

<br>

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

<br><br>
#### 2.2. *main.nf* file in a public github / gitlab repository

Run the following command from where you want the results:

<pre>
nextflow run -hub pasteur gmillot/repertoire_profiler -r v1.0.0
</pre>

<br><br>

### 3. Distant running (example with the Pasteur cluster)

#### 3.1. Pre-execution

Go into the directory where the main.nf and nextflow.config files are.
Copy-paste this code:

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

<br><br>

#### 3.2. *main.nf* file in a cluster folder

Modify the second line of the code below, and run from where the *main.nf* and *nextflow.config* files are (which has been set thanks to the EXEC_PATH variable above):

<pre>
HOME_INI=$HOME
HOME="${HELIXHOME}/repertoire_profiler/" # $HOME changed to allow the creation of .nextflow into /$HELIXHOME/repertoire_profiler/, for instance. See NFX_HOME in the nextflow software script
nextflow run main.nf -c nextflow.config # or nextflow run main.nf -c nextflow.config --modules ${MODULES} in order to have all the used module versions recorded into the report file 
HOME=$HOME_INI
</pre>

<br><br>

#### 3.3. *main.nf* file in the public gitlab repository

Modify the first and third lines of the code below, and run (results will be where the EXEC_PATH variable has been set above):

<pre>
VERSION="v1.0"
HOME_INI=$HOME
HOME="${HELIXHOME}/repertoire_profiler/" # $HOME changed to allow the creation of .nextflow into /$HELIXHOME/repertoire_profiler/, for instance. See NFX_HOME in the nextflow software script
nextflow run -hub pasteur gmillot/repertoire_profiler -r $VERSION -c $HOME/nextflow.config
HOME=$HOME_INI
</pre>

<br><br>
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

Register at Itol as explained in the [Prerequisite](#1-prerequisite) section above, or set the `phylo_tree_itol_subscription` parameter of the *nextflow.config* file to `FALSE` and rerun.

<br><br>
## OUTPUT

By default, all the results are returned in a *result* folder where the *main.nf* executed file is located (created if does not exist). This can be changed using the `out_path_ini` parameter of the *nextflow.config* file. By default, each execution produces a new folder named *repertoire_profiler_\<ID\>*, created inside the *result* folder and containing all the outputs of the execution. The name of the folder can be changed using the `result_folder_name` parameter of the *nextflow.config* file. The new name file will be followed by an \<ID\> in all cases.
<br><br>
An example of results obtained with the dataset is present at this address: https://zenodo.org/records/15132203/files/repertoire_profiler_1743690584.zip.
<br><br>
Complete informations are in the Protocol 144-rev0 Ig clustering - Immcantation.docx (contact Gael Millot).
<br><br>
If the text is cut in the table, reload the page or change the width of the window.
<br><br>
<div style="overflow-x:auto; max-width:100%;">
<table style="width:100%; border-collapse:collapse; overflow-wrap: anywhere;table-layout:fixed; word-break:break-all;">
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            repertoire_profiler_&lt;UNIQUE_ID&gt; folders and files
        </th>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            Description
        </th>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            report.html
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            HTML report presenting the main results of the data processing. Its content will be gradually expanded to provide a more comprehensive analysis. 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            reports
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            Folder containing all the reports of the different processes, as well as the <i>nextflow.config</i> file used and the map of the processes execution (<i>nf_dag.png</i> file). 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            repertoires
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            Folder containing the repertoires, i.e., contingency tables of the V, J and C (constant) allele and gene usage from the <i>productive_seq.tsv</i> file (see below). Warning: the script currently takes the first annotation of the imgt annotation if several are presents in the v_call, j_call and c_call columns of the <i>productive_seq.tsv</i> file. (e.g., v_call with IGKV1-39*01,IGKV1D-39*01), so that contingencies are identical to those from the donut frequencies, that use germline_v_call and germline_j_call columns (allele reassignment by the <code>CreateGermlines.py</code> tool of immcantation) 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            figures
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            Folder containing:<ul style="padding-left:1.2em; margin:0;"><li><b>png</b>: folder containing the graphs in png format. </li><li><b>svg</b>: folder containing the graphs in svg vectorial format.</ul> 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            RData
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            Folder containing, for each clonal group, objects that can be used in R to further analyze or plot the data:
            <br><ul style="padding-left:1.2em; margin:0;"><li>db: tibble data frame resulting from the import by the <code>alakazam::readChangeoDb()</code> function.
            <br></li><li>clones: db in the airClone format.
            <br></li><li>trees: output of the <code>dowser::getTrees()</code> function using the clones object as input (igphylm tree).
            </ul><br>Also contains the all_trees.RData file that combine the trees R objects of the different files in a single trees object. 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            alignments
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-word;
    max-width:100%; overflow-wrap:anywhere;">
            Folder containing alignment files in amino-acid (aa folder) and nucleotide (nuc folder) sequences: 
            <ul style="padding-left:1.2em; margin:0;word-break:break-word;
    max-width:100%;"><li>*.html: visualization of the alignments. For clonal groups, each file is named as <i>&lt;KIND&gt;&#8203;_clone_id_&#8203;&lt;CLONE_ID&gt;&#8203;_&lt;V_GENE&gt;&#8203;_&lt;J_GENE&gt;&#8203;_aligned_&lt;&#8203;nuc OR aa&gt;&#8203;.html</i>. For all sequences, the name is <i>&lt;KIND&gt;&#8203;_aligned_&lt;&#8203;nuc OR aa&gt;&#8203;.html</i>.
            </li><li>*.fasta: aligned sequences. Each file is named as the corresponding <i>.html</i> file.
            </li><li>*.gff: file used to add domain features in the corresponding <i>.html</i> file, and named as this one.</ul><br>
            See the tsv/productive_seq.tsv file for the sequences used in all aligned sequences, and the tsv/clone_assigned_seq.tsv file for the repartition of the sequences by clonal ID (clone_id column).<br><br>
            Warning: alignments are perfomed using <a href="https://mafft.cbrc.jp/alignment/server/index.html">Mafft</a> (no constraint on the V, D, J regions).<br><br>
            Warning: html, fasta and gff files can be absent for clonal groups if none have at least the number of sequences indicated by the `clone_nb_seq` parameter in the <i>nextflow.config</i> file.
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            phylo
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            Folder containing phylogenic tree files in amino-acid (aa folder) and nucleotide (nuc folder) sequences of each .fasta aligment file present in the <i>alignments</i> folder. Trees are obtained using `iqtree -m phylo_tree_model_file+GTR+I+R6` for nucleotide sequences and `iqtree -m phylo_tree_model_file+I+R6` for aa sequences, with `phylo_tree_model_file` a parameter of the <i>nextflow.config</i> file. By default, this parameter is a <a href="https://doi.org/10.1093/molbev/msu340">AB model</a> file dedicated to antibodies.
            <ul style="padding-left:1.2em; margin:0;"><li>*fasta.treefile: phylogenetic tree in Newick format for each sequence group in the phylo file (not constructed for groups with fewer than 4 sequences).</li><li>*fasta_itol_url.txt (if the `phylo_tree_itol_subscription` parameter of the <i>newtflow.config</i> file is `TRUE`): URL of the phylogenetic tree in iTOL.</li></ul>
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            fasta
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            Folder containing all fasta files:
            <ul style="padding-left:1.2em; margin:0;word-break:break-all; overflow-wrap:anywhere;"><li><b>for_alignment_nuc</b>: folder containing the nucleotide sequences belonging to the same clonal group (i.e., with the same value in the <i>clone_id</i> column of the  <i>clone_assigned_seq.tsv</i> file).<br>Warning: fasta files can be absent for clonal groups if none have at least the number of sequences indicated by the `clone_nb_seq` parameter in the <i>nextflow.config</i> file.<br>The kind of sequences depends on what has been selected in the `align_kind` parameter of the <i>nextflow.config</i> file. Each file is named as <i>&lt;KIND&gt;&#8203;_clone_id_&lt;&#8203;CLONE_ID&gt;_&lt;&#8203;V_GENE&gt;_&lt;&#8203;J_GENE&gt;&#8203;.fasta</i>.<br>If <i>align_kind</i> is "query", "igblast_full" or "trimmed", then the sequence of the <i>germline_alignment</i> column is also added with the following rule: 1) if all the <i>germline_alignment</i> sequences are not identical in the clonal group, the most frequent one is used and 2) if there is more than one most frequent germline sequence, then the first one is taken. See the <i>Tsv2fastaGff.log</i> file in the <i>report</i> folder. Warning: the germline sequence is always the variable region only (thus, can be shorter compared to the other sequences).<br>If <i>align_kind</i> is "sequence_alignment", "v_sequence_alignment", "d_sequence_alignment", "j_sequence_alignment" or "c_sequence_alignment", then the sequence of the corresponding <i>germline_alignment</i> column is also added with the same rules as above. Sequences are made of ATGCN.
            </li><li><b>for_alignment_aa</b>: idem to <i>for_alignment_nuc</i> but using the corresponding aa columns of the <i>clone_assigned_seq.tsv</i> file. Each file is named as <i>&lt;KIND&gt;&#8203;_aa_clone_id_&lt;&#8203;CLONE_ID&gt;_&lt;&#8203;V_GENE&gt;_&lt;&#8203;J_GENE&gt;&#8203;.fasta</i>. Germline sequence addition follows the same rules as in *for_alignment_nuc*, but using the corresponding aa column.Thus, evan if the name of the germline sequence is the same in both the nuc and aa alignment fasta file (for instance, germline&#8203;_alignment&#8203;_clone&#8203;_id&#8203;_10&#8203;_IGHV4-34&#8203;_IGHJ6), while sequences comes from 2 different columns (for instance, *germline_alignment* for nuc and *germline_alignment_aa* for aa). Sequences are made of single letter aa with X and * added.
            </li><li><b>productive_nuc</b>: folder containing fasta files of the nucleotidic sequences in the <i>productive_seq.tsv</i> file.
            <ul style="padding-left:1.2em; margin:0;"><li><b>query</b>: folder of the sequences from the <i>sequence</i> column. Sequences are made of ATGCN.
            </li><li><b>removed</b>: folder of the sequences from the <i>removed_sequence</i> column. Sequences are made of ATGCN.
            </li><li><b>trimmed</b>: folder of the sequences from the <i>trimmed_sequence</i> column. Sequences are made of ATGCN. Warning: empty fasta means that trimming was not possible by <code>seqkit subseq</code> because sequences in the <i>sequence_alignment</i> and <i>sequence_alignment_with_gaps</i> columns are not identical (align_seq_identical column is FALSE).
            </li><li><b>align</b>: folder of the sequences from the <i>sequence_alignment</i> column. Sequences are made of ATGCN-.
            </li><li><b>align_with_gaps</b>: folder of the sequences from the <i>sequence_alignment_with_gaps</i> column. Sequences are made of ATGCN.-.
            </li></ul></li><li><b>productive_aa</b>: folder containing fasta files of the translated sequences in the <i>productive_seq.tsv</i> file.
            <ul style="padding-left:1.2em; margin:0;"><li><b>query</b>: folder of the sequences from the <i>query_sequence_aa</i> column. Sequences are made of single letter aa with X and * added.
            </li><li><b>igblast</b>: folder of the sequences from the <i>sequence_aa</i> column. Sequences are made of single letter aa with X and * added.
            </li><li><b>trimmed</b>: folder of the sequences from the <i>trimmed_sequence_aa</i> column. Sequences are made of single letter aa with X and * added. Warning: empty fasta means that trimming was not possible by <code>seqkit subseq</code> because sequences in the <i>sequence_alignment</i> and <i>sequence_alignment_with_gaps</i> columns are not identical (align_seq_identical column is FALSE).
            </li><li><b>align</b>: folder of the sequences from the <i>sequence_alignment_aa</i> column. Sequences are made of single letter aa with X, * and - added.</li></ul>
            <br><b>IMPORTANT NOTE ABOUT CLONAL GROUPS</b> :<br>If the `clone_strategy` parameter in <i>nextflow.config</i> is set to "set", a single clonal group may contain different gene assignments. In such cases, the germline sequence included in the aligned FASTA file is selected as the most frequent germline sequence within the group.<br>If multiple germline sequences are equally frequent, the one associated with the first sequence (after sorting the group by sequence ID) is used. A warning is written to the log file whenever multiple germline sequences are detected in a group. 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            pdf
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            Folder containing all the  <b>pdf files</b> described below : 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            - seq_distance.pdf
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            Distribution of the distances between the two nearest sequences (see the <i>nearest_distance</i> column in the <i>productive_seq.tsv</i> file). 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            - donuts.pdf
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            Donut plots showing the frequency of sequences among:<br><ul style="padding-left:1.2em; margin:0;"><li><b>all</b>: all the productive sequences (<i>productive_seq.tsv</i> output file).
            <br></li><li><b>annotated</b>: as the "all" donut but using all the productive sequences that have been annotated using the `meta_name_replacement` parameter of the <i>nextflow.config</i> file if not "NULL".
            <br></li><li><b>trees</b>: all the sequences used for germline trees (<i>germ_tree_seq.tsv</i> output file).</li>
            </ul>
            <br>Warning: the script currently takes the first annotation of the imgt annotation if several are presents in the v_call, j_call or c_call columns of the <i>productive_seq.tsv</i> and <i>germ_tree_seq.tsv</i> files.
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            - *_repertoire.pdf
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            heatmap of the files from the <i>repertoires</i> folder (see above), showing the frequency of 1) alleles or 2) genes used among all the all productive sequences ("all"), non empty cells ("non-zero") and "annotated" sequences (if metadata are provided). Non-zero means that unused alleles are removed from the heatmap (empty row or column). Warning: to build the repertoire contingencies, the script currently takes the first annotation of the imgt annotation if several are presents in the v_call, j_call or c_call columns of the <i>productive_seq.tsv</i> file. 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            tsv
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            Folder containing all the  <b>tsv files</b> described below : 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            - passed_igblast_seq.tsv
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            sequences annotated by igblast (more precisely by <code>AssignGenes.py igblast --format airr</code>). If empty, generate a subsequent nextflow failure. The number of lines in <i>failed_igblast_seq.tsv</i> and <i>passed_igblast_seq.tsv</i> is equal to the number of submitted .fasta files. 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            - failed_igblast_seq.tsv
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            sequences that failed to be annotated by igblast (empty file if all the sequences are annotated). 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            - productive_seq.tsv
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            Productive sequences. Productive <a href="https://docs.airr-community.org/en/stable/datarep/rearrangements.html#productive">means</a>: (1) coding region has an open reading frame, (2) no defect in the start codon, splicing sites or regulatory elements, (3) no internal stop codons, (4) an in-frame <a href="https://docs.airr-community.org/en/stable/datarep/rearrangements.html#junction-versus-cdr3">junction</a> region. See the <i>unproductive_seq.tsv</i> file for sequences that failed to be productive.<br><b>Names in bold</b> are columns that have been added or modified, compared to the output of <code>AssignGenes.py igblast --format airr</code>.<br><a href="https://docs.airr-community.org/en/stable/datarep/rearrangements.html#coordinate-numbering">Coordinates</a> are 1-based numbering with closed intervals. This means that 1) the very first position in sequence is numbered 1 (not 0) and 2) both the start and the end positions are included in the range.<br>Column description (from <a href="https://docs.airr-community.org/en/stable/datarep/rearrangements.html#fields">here</a> for the names not in bold, with fields not kept <a href="./bin/fields_not_kept.txt">here</a>):
            <br><ul style="padding-left:1.2em; margin:0;"><li><b>sequence_id</b>: equivalent to <i>initial_sequence_id</i> column but modified with new names according to the `meta_name_replacement` parameter of the <i>nextflow.config</i> file.
            <br></li><li><b>initial_sequence_id</b>: optional. Only present if the <i>meta_path</i> and `meta_name_replacement` parameters of the <i>nextflow.config</i> file are non NULL. Originaly the <i>sequence_id</i> column of the <code>AssignGenes.py igblast --format airr</code> output. Unique query sequence identifier for the Rearrangement. Most often this will be the input sequence header or a substring thereof, but may also be a custom identifier defined by the tool in cases where query sequences have been combined in some fashion prior to alignment. When downloaded from an AIRR Data Commons repository, this will usually be a universally unique record locator for linking with other objects in the AIRR Data Model.
            <br></li><li>sequence: the query (input) nucleotide sequence. Usually, this is the unmodified input sequence, which may be reverse complemented if necessary. In some cases, this field may contain consensus sequences or other types of collapsed input sequences if these steps are performed prior to alignment.
            <br></li><li>sequence_aa: Translation in aa of the query nucleotide sequence performed by <code>AssignGenes.py igblast --format airr</code>. Warning: see the <i>aa_identical</i> column description.
            <br></li><li>stop_codon: True if the aligned sequence contains a stop codon.
            <br></li><li>vj_in_frame: True if the V and J gene alignments are in-frame.
            <br></li><li>v_frameshift: True if the V gene in the query nucleotide sequence contains a translational <a href="https://docs.airr-community.org/en/stable/datarep/rearrangements.html#frameshifts">frameshift</a> relative to the frame of the V gene reference sequence.
            <br></li><li>productive: True if the V(D)J sequence is predicted to be productive.
            <br></li><li>rev_comp: True if the alignment is on the opposite strand (reverse complemented) with respect to the query sequence. If True then all output data, such as alignment coordinates and sequences, are based on the reverse complement of 'sequence'.
            <br></li><li>complete_vdj: True if the sequence alignment spans the entire V(D)J region. Meaning, sequence_alignment includes both the first V gene codon that encodes the mature polypeptide chain (i.e., after the leader sequence) and the last complete codon of the J gene (i.e., before the J-C splice site). This does not require an absence of deletions within the internal FWR and CDR regions of the alignment.
            <br></li><li>d_frame: Numerical reading frame (1, 2, 3) of the D gene in the query nucleotide sequence, where frame 1 is relative to the first codon of D gene reference sequence.
            <br></li><li>v_call: V gene (i.e., cassette) with allele after the star (e.g., IGHV4-59*01 if using IMGT/GENE-DB). Sometimes, Igblast cannot distinguished between several reference sequences (cassette or allele), and are all provided, comma separated (e.g., IGHJ5*01,IGHJ5*02).d_call: as in the <i>v_call</i> column but for D gene with allele.
            <br></li><li>j_call: as in the <i>v_call</i> column but for J gene with allele.
            <br></li><li>c_call: as in the <i>v_call</i> column but for Constant region gene with allele.
            <br></li><li>sequence_alignment: Aligned portion of query sequence (i.e., <i>sequence</i> column). Typically, this will include  <b>only the V(D)J region</b>, but that is not a requirement, explaining why this sequence lacks the 3' part of the sequence indicated in the <i>sequence</i> column. Warning, it seems that <code>AssignGenes.py igblast --format airr</code> (currently used in this pipeline) does not introduce the IMGT-gaps spacers as with <code>AssignGenes.py igblast --format blast</code> (previously used in this pipeline). See the <i>sequence_alignment_with_gaps</i> column to have the IMGT-gaps spacers (but see also the <i>align_seq_identical</i> column to check that both sequences are identical, gaps excluded). Of note, <i>v_sequence_alignment</i>, <i>d_sequence_alignment</i>, <i>j_sequence_alignment</i> and <i>c_sequence_alignment</i> sequence assembly is not always equivalent to <i>sequence_alignment</i>, due to extra nucleotides between v-d, d-j and j-c, and potential overlap between j-c. Sequences are made of ATGCN-.
            <br></li><li>germline_alignment: Assembled, aligned, full-length inferred germline sequence spanning the same region as the sequence_alignment field (typically the V(D)J region) and including the same set of corrections and IMGT-gaps spacers (if any). Thus, If well understood, this sequence is built from V, D, J reference sequences from the IMGT database that match the sequence in sequence_alignment, with potential IMGT-gaps spacers already present in the sequences of the database kept. The aligned parts of these V, D, J references are stitched together, and the <code>N</code> nucleotides in the germline alignment sequence are positions that belong to no templated nucleotides, added during V(D)J recombination by enzymes involved in this process, and not part of the reference. These 'N's only appear in the junction region where DJ recombination and V-DJ recombination happened and are present where VDJ reference cassettes don't match each other. Warning: sequences with the same clone_ID can have different germline_aligments since the column does not come from <code>CreateGermlines.py</code>. Warning, it seems that <code>AssignGenes.py igblast --format airr</code> (currently used in this pipeline) does not introduce the IMGT-gaps spacers as with <code>AssignGenes.py igblast --format blast</code> (previously used in this pipeline). See the <i>germline_alignment_with_gaps</i> column to have the IMGT-gaps spacers. Of note, <i>v_germline_alignment</i>, <i>d_germline_alignment</i>, <i>j_germline_alignment</i> and <i>c_germline_alignment</i> sequence assembly is not always equivalent to <i>germline_alignment</i>, due to extra nucleotides between v-d, d-j and j-c, and potential overlap between j-c. Sequences are made of ATGCN-.
            <br></li><li>sequence_alignment_aa: Amino acid translation of the <i>sequence_alignment</i> column.  Warning: see the <i>aa_identical</i> column description. Of note, <i>v_sequence_alignment_aa</i>, <i>d_sequence_alignment_aa</i>, <i>j_sequence_alignment_aa</i> and <i>c_sequence_alignment_aa</i> sequence assembly is not always equivalent to <i>sequence_alignment_aa</i>, due to extra amino acids between v-d, d-j and j-c. Sequences are made of single letter aa with X, * and - added.
            <br></li><li>germline_alignment_aa: Amino acid translation of the <i>germline_alignment</i> column. Warning: see the <i>aa_identical</i> column description. Of note, <i>v_germline_alignment_aa</i>, <i>d_germline_alignment_aa</i>, <i>j_germline_alignment_aa</i> and <i>c_germline_alignment_aa</i> sequence assembly is not always equivalent to <i>germline_alignment_aa</i>, due to extra amino acids between v-d, d-j and j-c. Sequences are made of single letter aa with X, * and - added.
            <br></li><li>v_alignment_start: Start position of the V gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).
            <br></li><li>v_alignment_end: End position of the V gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).
            <br></li><li>d_alignment_start: Start position of the D gene in both the sequence_alignment and germline_alignment fields (1-based closed interval). Warning: D can overlap V.
            <br></li><li>d_alignment_end: End position of the D gene in both the sequence_alignment and germline_alignment fields (1-based closed interval). Warning: D can overlap J.
            <br></li><li>j_alignment_start: Start position of the J gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).
            <br></li><li>j_alignment_end: End position of the J gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).
            <br></li><li>c_alignment_start: Start position of the C gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).
            <br></li><li>c_alignment_end: End position of the C gene alignment in both the sequence_alignment and germline_alignment fields (1-based closed interval).
            <br></li><li>v_sequence_alignment: Aligned portion of query sequence (i.e., <i>sequence</i> column) assigned to the V gene, including any indel corrections or numbering spacers. Sequences are made of ATGC-.
            <br></li><li>v_sequence_alignment_aa: Amino acid translation of the v_sequence_alignment field. Sequences are made of single letter aa with - added.
            <br></li><li>v_germline_alignment: Aligned V gene germline sequence spanning the same region as the v_sequence_alignment field and including the same set of corrections and spacers (if any). Sequences are made of ATGC-.
            <br></li><li>v_germline_alignment_aa: Amino acid translation of the v_germline_alignment field.
            <br></li><li>d_sequence_alignment: Aligned portion of query sequence (i.e., <i>sequence</i> column) assigned to the D gene, including any indel corrections or numbering spacers. Sequences are made of ATGC-.
            <br></li><li>d_sequence_alignment_aa: Amino acid translation of the d_sequence_alignment field.
            <br></li><li>d_germline_alignment: Aligned D gene germline sequence spanning the same region as the d_sequence_alignment field and including the same set of corrections and spacers (if any). Sequences are made of single letter aa with - added.
            <br></li><li>d_germline_alignment_aa: Amino acid translation of the d_germline_alignment field.
            <br></li><li>j_sequence_alignment: Aligned portion of query sequence (i.e., <i>sequence</i> column) assigned to the J gene, including any indel corrections or numbering spacers. Sequences are made of ATGC-.
            <br></li><li>j_sequence_alignment_aa: Amino acid translation of the j_sequence_alignment field.
            <br></li><li>j_germline_alignment: Aligned J gene germline sequence spanning the same region as the j_sequence_alignment field and including the same set of corrections and spacers (if any). Sequences are made of single letter aa with - added.
            <br></li><li>j_germline_alignment_aa: Amino acid translation of the j_germline_alignment field.
            <br></li><li>c_sequence_alignment: Aligned portion of query sequence (i.e., <i>sequence</i> column) assigned to the constant region, including any indel corrections or numbering spacers. Sequences are made of ATGC-.
            <br></li><li>c_sequence_alignment_aa: Amino acid translation of the c_sequence_alignment field. Sequences are made of single letter aa with - added.
            <br></li><li>c_germline_alignment: Aligned constant region germline sequence spanning the same region as the c_sequence_alignment field and including the same set of corrections and spacers (if any).
            <br></li><li>c_germline_alignment_aa: Amino acid translation of the c_germline_aligment field. Sequences are made of single letter aa with - added.
            <br></li><li>fwr1: Nucleotide sequence of the aligned FWR1 region of the query sequence (i.e., FWR1 region of the <i>sequence_alignment</i> field, which is the input sequence with IMGT gaps added). Sequences are made of ATGC.
            <br></li><li>fwr1_aa: Amino acid translation of the fwr1 field. Sequences are made of single letter aa.
            <br></li><li>cdr1: Nucleotide sequence of the aligned CDR1 region of the query sequence (i.e., CDR1 region of the <i>sequence_alignment</i> field, which is the input sequence with IMGT gaps added). Sequences are made of ATGC.
            <br></li><li>cdr1_aa: Amino acid translation of the cdr1 field. Sequences are made of single letter aa.
            <br></li><li>fwr2: Nucleotide sequence of the aligned FWR2 region of the query sequence (i.e., FWR2 region of the <i>sequence_alignment</i> field, which is the input sequence with IMGT gaps added). Sequences are made of ATGC.
            <br></li><li>fwr2_aa: Amino acid translation of the fwr2 field. Sequences are made of single letter aa.
            <br></li><li>cdr2: Nucleotide sequence of the aligned CDR2 region of the query sequence (i.e., CDR2 region of the <i>sequence_alignment</i> field, which is the input sequence with IMGT gaps added). Sequences are made of ATGC.
            <br></li><li>cdr2_aa: Amino acid translation of the cdr2 field. Sequences are made of single letter aa.
            <br></li><li>fwr3: Nucleotide sequence of the aligned FWR3 region of the query sequence (i.e., FWR3 region of the <i>sequence_alignment</i> field, which is the input sequence with IMGT gaps added). Sequences are made of ATGC.
            <br></li><li>fwr3_aa: Amino acid translation of the fwr3 field. Sequences are made of single letter aa.
            <br></li><li>fwr4: Nucleotide sequence of the aligned FWR4 region of the query sequence (i.e., FWR4 region of the <i>sequence_alignment</i> field, which is the input sequence with IMGT gaps added). Sequences are made of ATGC.
            <br></li><li>fwr4_aa: Amino acid translation of the fwr4 field. Sequences are made of single letter aa.
            <br></li><li>cdr3: Nucleotide sequence of the aligned CDR3 region of the query sequence (i.e., CDR3 region of the <i>sequence_alignment</i> field, which is the input sequence with IMGT gaps added). Sequences are made of ATGC.
            <br></li><li>cdr3_aa: Amino acid translation of the cdr3 field. Sequences are made of single letter aa.
            <br></li><li>junction: Junction region nucleotide sequence, where the <a href="https://docs.airr-community.org/en/stable/datarep/rearrangements.html#junction-versus-cdr3">Junction</a> is defined as the CDR3 plus the two flanking conserved codons (including the conserved cysteine and tryptophan/phenylalanine residues, while CDR3 excludes those). Sequences are made of ATGC.
            <br></li><li>junction_length: Number of nucleotides in the junction sequence.
            <br></li><li>junction_aa: Amino acid translation of the junction. Sequences are made of single letter aa.
            <br></li><li>junction_aa_length: Number of amino acids in the junction sequence.
            <br></li><li>v_score: Alignment score for the V gene.
            <br></li><li>d_score: Alignment score for the D gene alignment.
            <br></li><li>j_score: Alignment score for the J gene alignment.
            <br></li><li>c_score: Alignment score for the C gene alignment.
            <br></li><li>v_cigar: CIGAR string for the V gene alignment.
            <br></li><li>d_cigar: CIGAR string for the D gene alignment.
            <br></li><li>j_cigar: CIGAR string for the J gene alignment.
            <br></li><li>c_cigar: CIGAR string for the C gene alignment.
            <br></li><li>v_support: V gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the V gene assignment as defined by the alignment tool.
            <br></li><li>d_support: D gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the D gene as defined by the alignment tool.
            <br></li><li>j_support: J gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the J gene assignment as defined by the alignment tool.
            <br></li><li>c_support: C gene alignment E-value, p-value, likelihood, probability or other similar measure of support for the C gene assignment as defined by the alignment tool.
            <br></li><li>v_identity: Fractional identity for the V gene alignment (proportion).
            <br></li><li>d_identity: Fractional identity for the D gene alignment (proportion).
            <br></li><li>j_identity: Fractional identity for the J gene alignment (proportion).
            <br></li><li>c_identity: Fractional identity for the C gene alignment (proportion).
            <br></li><li>v_sequence_start: Start position of the V gene in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>v_sequence_end: End position of the V gene in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>v_germline_start: Alignment start position in the V gene reference sequence (1-based closed interval).
            <br></li><li>v_germline_end: Alignment end position in the V gene reference sequence (1-based closed interval).
            <br></li><li>d_sequence_start: Start position of the D gene in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>d_sequence_end: End position of the D gene in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>d_germline_start: Alignment start position in the D gene reference sequence for the D gene (1-based closed interval).
            <br></li><li>d_germline_end: Alignment end position in the D gene reference sequence for the D gene (1-based closed interval).
            <br></li><li>j_sequence_start: Start position of the J gene in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>j_sequence_end: End position of the J gene in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>j_germline_start: Alignment start position in the J gene reference sequence (1-based closed interval).
            <br></li><li>j_germline_end: Alignment end position in the J gene reference sequence (1-based closed interval).
            <br></li><li>c_sequence_start: Start position of the C gene in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>c_sequence_end: End position of the C gene in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>c_germline_start: Alignment start position in the C gene reference sequence (1-based closed interval).
            <br></li><li>c_germline_end: Alignment end position in the C gene reference sequence (1-based closed interval).
            <br></li><li>fwr1_start: FWR1 start position in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>fwr1_end: FWR1 end position in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>cdr1_start: CDR1 start position in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>cdr1_end: CDR1 end position in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>fwr2_start: FWR2 start position in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>fwr2_end: FWR2 end position in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>cdr2_start: CDR2 start position in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>cdr2_end: CDR2 end position in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>fwr3_start: FWR3 start position in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>fwr3_end: FWR3 end position in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>fwr4_start: FWR4 start position in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>fwr4_end: FWR4 end position in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>cdr3_start: CDR3 start position in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>cdr3_end: CDR3 end position in the query sequence (i.e., <i>sequence</i> column, 1-based closed interval).
            <br></li><li>np1: Nucleotide sequence of the combined N/P region between the V gene and first D gene alignment or between the V gene and J gene alignments.
            <br></li><li>np1_length: Number of nucleotides between the V gene and first D gene alignments or between the V gene and J gene alignments.
            <br></li><li>np2: Nucleotide sequence of the combined N/P region between either the first D gene and J gene alignments or the first D gene and second D gene alignments.
            <br></li><li>np2_length: Number of nucleotides between either the first D gene and J gene alignments or the first D gene and second D gene alignments.
            <br></li><li><b>sequence_alignment_with_gaps</b>: Aligned portion of query sequence (i.e., <i>sequence</i> column), including any indel corrections or numbering spacers, such as IMGT-gaps (explained <a href="https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html">here</a>), in order to fit the <a href="https://www.imgt.org/3Dstructure-DB/doc/IMGTCollier-de-Perles.shtml">collier de perle represention</a>. Typically, this will include  <b>only the V(D)J region</b>, but that is not a requirement, explaining why this sequence lacks the 3' part of the sequence indicated in the <i>sequence</i> column. Of note, these sequences correspond to the <i>sequence_alignment</i> column of the .tsv ouptut obtained with <code>AssignGenes.py igblast --format blast</code> and <code>MakeDb.py igblast</code>. Sequences are made of ATGC.-.
            <br></li><li><b>germline_alignment_with_gaps</b>: Assembled, aligned, full-length inferred germline sequence spanning the same region as the sequence_alignment field (typically the V(D)J region) and including the same set of corrections and IMGT-gaps spacers (if any). Thus, If well understood, this sequence is built from V, D, J reference sequences from the IMGT database that match the sequence in sequence_alignment, with potential IMGT-gaps spacers already present in the sequences of the database kept. The aligned parts of these V, D, J references are stitched together, and the <code>N</code> nucleotides in the germline alignment sequence are positions that belong to no templated nucleotides, added during V(D)J recombination by enzymes involved in this process, and not part of the reference. These 'N's only appear in the junction region where DJ recombination and V-DJ recombination happened and are present where VDJ reference cassettes don't match each other. Warning: sequences with the same clone_ID can have different germline_aligments since the column does not come from <code>CreateGermlines.py</code>. Of note, correspond to the <i>sequence_alignment</i> column obtained with <code>AssignGenes.py igblast --format blast</code>. Sequences are made of ATGCN.-.
            <br></li><li><b>trimmed_sequence</b>: identical to the query sequence (in the <i>sequence</i> column), except a potential trimming of the 5' part that removes all the nucleotides before the fwr1 region (corresponding to the leader peptide). This should guarantees a translation of these sequences without premature stops. Trimming is performed by comparing the 5' part of the <i>sequence_alignment</i> and the <i>sequence</i> sequences. Sequences are made of ATGCN.
            <br></li><li><b>is_sequence_trimmed</b>: indicates if the sequences in the <i>sequence</i> column are trimmed or not. If FALSE, then the data in both <i>sequence</i> and <i>sequence</i> should be identical (<i>removed_sequence</i> column is NA), or that a problem of trimming has occured (<i>removed_sequence</i> column is not NA, see the <i>trimmed</i> column).
            <br></li><li><b>removed_sequence</b>: removed sequence in 5'. NA if no removal (<i>is_sequence_trimmed</i> column FALSE). Sequences are made of ATGCN.
            <br></li><li><b>trimmed_sequence_aa</b>: Translation in aa of the <i>sequence</i> column using <code>seqkit translate</code>. Sequences are made of single letter aa with X and * added.
            <br></li><li><b>query_sequence_aa</b>: Translation in aa of the <i>sequence</i> column (i.e., query nucleotide sequence) using <code>seqkit translate</code>. If <i>is_sequence_trimmed</i> column is TRUE, then some * are expected to be present in the sequence (stop codons), since the query sequence does not necessarily start with a start codon (primer sequence addition for instance). Sequences are made of single letter aa with X and * added.
            <br></li><li><b>align_seq_identical</b>: if TRUE, it means that both the sequences in the <i>sequence_alignment</i> (coming from <code>AssignGenes.py igblast --format airr</code>) and <i>sequence_alignment_with_gaps</i> (coming from <code>AssignGenes.py igblast --format blast</code> and <code>MakeDb.py igblast</code>) columns are identical, gaps (dots) excluded. If FALSE, it raises concerns that remain to be explained.
            <br></li><li><b>aa_identical</b>: if TRUE, it means that both the sequences in the <i>sequence_aa</i> and <i>trimmed_sequence_aa</i> columns are identical. If FALSE, it can be because of an extra aa at the end of the <i>sequence_aa</i> sequence, due to incomplete last codon but with aa inferred by <code>AssignGenes.py igblast --format airr</code>. Otherwise, it raises problems in the output of <code>AssignGenes.py igblast --format airr</code>, notably the correspondance between: 1) <i>sequence_alignment_aa</i> and <i>sequence_alignment</i> columns, 2) <i>germline_alignment_aa</i> and <i>germline_alignment</i> columns. In other words, translation of <i>sequence_alignment</i> and <i>germline_alignment</i> might not results in the corresponding columns.
            <br></li><li><b>sequence_aa_stop</b>: Is there any stop codon inside the <i>sequence_aa</i> sequence? TRUE if yes, FALSE if no. Warning: TRUE is not expected for productive sequences.
            <br></li><li><b>sequence_alignment_aa_stop</b>: Is there any stop codon inside the <i>sequence_alignment_aa</i> sequence? TRUE if yes, FALSE if no. Warning: TRUE is not expected for productive sequences.
            <br></li><li><b>germline_alignment_aa_stop</b>: Is there any stop codon inside the <i>germline_alignment_aa</i> sequence? TRUE if yes, FALSE if no. Warning: TRUE is not expected for productive sequences.
            <br></li><li><b>trimmed_sequence_aa_stop</b>: Is there any stop codon inside the <i>trimmed_sequence_aa</i> sequence? TRUE if yes, FALSE if no. If TRUE and stars are present in the end of the aa sequence, it might come from a bad sequencing quality.
            <br></li><li><b>query_sequence_aa_stop</b>: Is there any stop codon inside the <i>query_sequence_aa</i> sequence? TRUE if yes, FALSE if no.
            <br></li><li><b>&lt;OPTIONAL_COLUMN&gt;</b>: reporting the data from the `meta_legend` parameter of the <i>nextflow.config</i> file, if ever used. Example: "KD".
            <br></li><li><b>dist_nearest</b>: minimal distance from the nearest sequence using the model from the `clone_model` parameter (Haming by default). NA if no other sequences have same V, J and junction length or if another sequence is strictly identical (should be 0 but NA is returned). See <i>clone_id</i> below.
            <br></li><li><b>v_gene</b>: extracted from the <i>v_call</i> column but removing the allele specification after the *.
            <br></li><li><b>j_gene</b>: extracted from the <i>j_call</i> column but removing the allele specification after the *.
            <br></li><li><b>isotype_class</b>: extracted from the <i>c_call</i> column but indicating only the isotype.
            <br></li><li><b>c_gene</b>: extracted from the <i>c_call</i> column but removing the allele specification after the *.</li> 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            - unproductive_seq.tsv
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            Sequences that failed in productive annotations by igblast (empty file if all the sequences are productively annotated). 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            - clone_assigned_seq.tsv
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            Sequences from the <i>productive_seq.tsv</i> file with germline clustering (clone ID), allele reannotation (germinal_v_call and germinal_j_call columns) and mutation load added. Warning: the number of sequences (i.e., rows) can be lower than in the <i>productive_seq.tsv</i> file due to sequences that failed to be clone assigned (see the <i>non_clone_assigned_sequence.tsv</i> file).
            <br>Additional columns description (from <a href="https://docs.airr-community.org/en/stable/datarep/rearrangements.html#fields">here</a>):
            <br><ul style="padding-left:1.2em; margin:0;"><li>clone_id: Clone number. A same clone_id gathers all the sequences that putatively come from a same germline cell. See <a href="https://changeo.readthedocs.io/en/stable/examples/cloning.html#assigning-clones">here</a>, <a href="https://shazam.readthedocs.io/en/stable/vignettes/DistToNearest-Vignette/">here</a> and <a href="https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-015-0243-2">this article</a> for details. In summary: 1) grouping the sequences according to "same V, J and junction length" (to facilitate the distance computation), 2) for each group, 2x2 distance computation using by default (`clone_model` and `clone_normalize` parameters of the <i>nextflow.config</i> file) the <a href="https://biology.stackexchange.com/questions/23523/hamming-distance-between-two-dna-strings">Hamming distance</a>, 3) cutoff definition , 4) using the cutoff (`clone_distance` parameter of the <i>nextflow.config</i> file) to define clonal groups inside each "same V, J and junction length" group.
            <br></li><li>germline_alignment_d_mask: germline sequence of the clonal group reconstructed by `CreateGermlines.py` with D masked (i.e., replaced by N, in the middle of the CDR3). Because the D-segment call for B cell receptor alignments is often low confidence, the default germline format ("dmask" option of the `clone_germline_kind` parameters of the <i>nextflow.config</i> file) places Ns in the N/P and D-segments of the junction region rather than using the D-segment assigned during reference alignment. This can be modified to generate a complete germline ("full" option) or a V-segment only germline ("vonly" option). This sequence is the assembly of the sequences in the <i>germline_v_seq</i>, <i>germline_d_seq</i> and <i>germline_j_seq</i> but with overlap of D on V and J and with the mask of D. Warning: this germline sequence is not necessarily identical to the one in the <i>germline_alignment</i> column. Indeed, the first one is identical to all sequences of the same clonal group. The latter is the reconstruction of the IMGT calling of each sequence. 
            <br></li><li>germline_v_call: V germline cassette
            <br></li><li>germline_d_call: D germline cassette (usually NA)
            <br></li><li>germline_j_call: J germline cassette
            <br></li><li>germline_v_seq: nucleotide sequence of the V germline cassette
            <br></li><li>germline_v_seq_no_gaps: nucleotide sequence of the V germline cassette without IMGT gaps
            <br></li><li>germline_d_seq: nucleotide sequence of the D germline cassette
            <br></li><li>germline_d_seq_no_gaps: nucleotide sequence of the D germline cassette without IMGT gaps
            <br></li><li>germline_j_seq: nucleotide sequence of the J germline cassette
            <br></li><li>germline_j_seq_no_gaps: nucleotide sequence of the J germline cassette without IMGT gaps
            <br></li><li>germline_d_mask_no_gaps: <i>germline_alignment_d_mask</i> column without IMGT gaps
            <br></li><li>germline_d_mask_aa_no_gaps: translation of the <i>germline_d_mask_no_gaps</i> column into amino-acids
            <br></li><li>mu_count_*_r: number of replacement mutations in the region indicated by the `clone_mut_regionDefinition` parameter of the <i>nextflow.config</i> file. See details <a href="https://shazam.readthedocs.io/en/stable/topics/observedMutations/#value">here</a>).
            <br></li><li>mu_count_*_s: number of silent mutations.
            <br></li><li>mu_count: number of replacement and silent mutations (sum of the previous columns).
            <br></li><li>mu_freq_*_r: frequency of replacement mutations.
            <br></li><li>mu_freq_*_s: frequency of silent mutations.
            <br></li><li>mu_freq: frequency of replacement and silent mutations (sum of the previous columns).
            <br></li><li><b>germline_v_gene</b>: extracted from the <i>germline_v_call</i> column but removing the allele specification after the *.
            <br></li><li><b>germline_d_gene</b>: extracted from the <i>germline_d_call</i> column but removing the allele specification after the *.
            <br></li><li><b>germline_j_gene</b>: extracted from the <i>germline_j_call</i> column but removing the allele specification after the *.
            </li> 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            - non_clone_assigned_sequence.tsv
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            Productive sequences that failed to be assigned to a clone ID by <code>DefineClones.py</code> (empty file if all the sequences are assigned). See details  <a href="https://changeo.readthedocs.io/en/latest/methods/clustering.html">here</a> but failure reasons are not explained. 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            - clone_id_count.tsv
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            number of sequences in each clonal group in the <i>clone_assigned_seq.tsv</i> file. 
        </td>
    </tr>
    <tr>
        <th style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            - donut_stat.tsv
        </th>
        <td style="white-space:normal; text-align:left; word-break:break-all; overflow-wrap:anywhere;">
            stats associated to the <i>donuts.pdf</i> file. 
        </td>
    </tr>
</table>
</div>

<br><br>
## VERSIONS


The different releases are tagged [here](https://gitlab.pasteur.fr/gmillot/repertoire_profiler/-/tags).

<br><br>
## LICENCE


This package of scripts can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchandability or fitness for a particular purpose.
See the GNU General Public License for more details at https://www.gnu.org/licenses or in the Licence.txt attached file.


<br><br>
## CITATION


Version V10.3:

[Dejoux A, et al. Sci Transl Med. 2024](https://www.science.org/doi/10.1126/scitranslmed.ado4463)
<br>

[Dejoux A, et al. J Allergy Clin Immunol. 2025](https://www.jacionline.org/article/S0091-6749(25)00113-7/fulltext)

<br><br>
## CREDITS

[Pascal Chappert](https://www.institut-necker-enfants-malades.fr/index.php?menu=team&rubric=teamtabs&idfac=mahevas#chart), INSERM U1151 Institut Necker Enfants Malades, Paris, France

[Frédéric Lemoine](), Institut Pasteur, Université Paris Cité, Bioinformatics and Biostatistics Hub, 75015 Paris, France

[Chloé Taurel](), University of Technology of Compiègne, Compiègne, France

[Gael A. Millot](https://gitlab.pasteur.fr/gmillot), Institut Pasteur, Université Paris Cité, Bioinformatics and Biostatistics Hub, 75015 Paris, France

<br><br>
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

Special acknowledgement to the team of [Kenneth Hoehn](https://medicine.yale.edu/profile/kenneth-hoehn/), Yale School of Medicine, New Haven, CT, USA

<br><br>
## WHAT'S NEW IN

#### v21.5

- Apptainer image pulling with time limit increased.

#### v21.4

- Donut updated so that it uses saferDev.

#### v21.3

- Donut updated for better management of legends.

#### v21.2

- Now take into account that nucleotide alignment sequences can have --- inside.

#### v21.1

- New column added in productive_seq.tsv.

#### v21.0

- New dataset in the nextflow.config file. Columns added in productive_seq.tsv. Bugs fixed.

#### v20.2

- Bugs fixed.

#### v20.1

- Remain to finish the alignments and phylo. The rest (repertoires notably) are ok.

#### v20.0

- Now the alignments are fixed and are on the variable region only

#### v19.2

- Bug fixed

#### v19.1

- Bug fixed and README completed

#### v19.0

- Now sequences are trimmed in 5' for the leader peptide if exists, and two new columns in productive_seq.tsv sequence_ini is_sequence_trimmed"

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

seq_not_displayed.tsv file added to better understand the absence of seq in trees<br>
V and J alleles added in the donut legend


#### v4.0

tree_meta_path modified so that it can now be a 'NULL' path<br>
V and J alleles added in tree titles<br>
Duplicated sequences can be removed or not from trees


#### v3.5

Bug solved


#### v3.4

Two Donut charts added


#### v3.3

Empty channel solved<br>
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



