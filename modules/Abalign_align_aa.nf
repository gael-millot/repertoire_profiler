
// Align the amino-acidic sequences that are already in fasta files (grouped by clonal groups)
process Abalign_align_aa {
    label 'abalign'
    
    input:
    tuple path(fasta_nuc), path(fasta_aa), val(seq_kind) // parallelization expected (by clonal groups over align_clone_nb sequences)
    val igblast_organism
    val igblast_B_heavy_chain
    val align_abalign_options
    
    output:
    tuple path(fasta_nuc), path(fasta_aa), path("*_align_aa.fasta"), val(seq_kind), emit: aligned_aa_ch
    path "abalign_align_aa.log", emit: abalign_align_aa_log_ch
    
    script:
    if( ! (igblast_B_heavy_chain == "TRUE" || igblast_B_heavy_chain == "FALSE") ){
        error "\n\n========\n\nERROR IN Abalign_align_aa PROCESS\n\nINVALID heavy_chain PARAMETER:\n${heavy_chain}\nMUST BE EITHER \"TRUE\" OR \"FALSE\"\n\n========\n\n"
    }
    parms="-al"
    if(igblast_B_heavy_chain == "TRUE"){parms="-ah"}
    // Choose the species parameter for abalign
    switch (igblast_organism) {
        case "mouse":
            species = "MU"
            break
        case "human":
            species = "HS"
            break
        case "rabbit":
            species = "OC"
            break
        case "rat":
            species = "MM"
            break
        case "rhesus_monkey":
            species = "RM"
            break
        default:
            error "\n\n========\n\nERROR IN Abalign_align_aa PROCESS\n\nINVALID igblast_organism PARAMETER:\n${igblast_organism}\nMUST BE EITHER \"mouse\" OR \"human\" OR \"rabbit\" OR \"rat\" OR \"rhesus_monkey\"\n\n========\n\n"
    }
    """
    #!/bin/bash -ue
    set -o pipefail
    FILENAME=\$(basename -- ${fasta_aa}) # recover a file name without path
    echo -e "\\n\\n################################\\n\\n\$FILENAME\\n\\n################################\\n\\n" |& tee -a abalign_align_aa.log
    echo -e "WORKING FOLDER:\\n\$(pwd)\\n\\n" |& tee -a abalign_align_aa.log
    /bin/Abalign_V2_Linux_Term/Abalign ${align_abalign_options} -i ${fasta_aa} ${parms} ${fasta_aa.baseName}_align_aa.fasta -sp ${species}. |& tee -a abalign_align_aa.log || true
    # -g   (IMGT Numbering scheme, default)
    # -lc length.txt -lg 1 # single length computed spanning the indicated regions. For instance, -lc length.txt -lg 1,2,3,4,5,6,7 returns a single length
    """
}

