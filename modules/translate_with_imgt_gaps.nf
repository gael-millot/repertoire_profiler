
process Translate_with_IMGT_gaps {
    label 'immcantation'

    input:
    path add_dotted_coord_ch// parallelization 

    output:
    path "add_aa_imgt.tsv", emit: add_aa_imgt_ch
    path "Translate_with_IMGT_gaps.log", emit: translate_with_IMGT_gaps_ch_log

    script:
    """
python3 - <<'EOF'
import pandas as pd
from Bio.Seq import Seq

# --- define translation function (keeps dots and hyphens) ---
def translate_with_gaps(nuc_seq):
    protein_seq = []
    seq_len = len(nuc_seq)
    valid_bases = set("ATGCatgc")
    for i in range(0, seq_len, 3):
        codon = nuc_seq[i:i+3]
        # incomplete codon at end → gap placeholder
        if len(codon) < 3:
            protein_seq.append('-')
            continue

        # all dots
        if all(b == '.' for b in codon):
            protein_seq.append('.')
            continue

        # all hyphens
        if all(b == '-' for b in codon):
            protein_seq.append('-')
            continue

        # if any invalid base, dot, or hyphen inside codon → treat as hyphen
        if any(b not in valid_bases for b in codon):
            protein_seq.append('-')
        else:
            try:
                aa = str(Seq(codon.upper()).translate(to_stop=False))
            except Exception:
                aa = 'X'
            protein_seq.append(aa)
    return ''.join(protein_seq)

# --- 1) import TSV file (header + single line) ---
df = pd.read_csv("${add_dotted_coord_ch}", sep="\\t")

# --- 2) get nucleotide sequence and 3) translate keeping dots and hyphens ---
df["sequence_alignment_with_gaps_aa"] = df["sequence_alignment_with_gaps"].apply(translate_with_gaps)

# --- 5) compare AA sequences (ignore dots/hyphens) ---
def compare_ignore_gaps(seq1, seq2):
    s1 = seq1.replace('.', '').replace('-', '')
    s2 = seq2.replace('.', '').replace('-', '')
    return s1 == s2

df["sequence_alignment_aa_identical"] = df.apply(
    lambda x: compare_ignore_gaps(x["sequence_alignment_with_gaps_aa"], x["sequence_alignment_aa"]),
    axis=1
)

# --- 6 & 7) insert new columns in specific positions ---
# insert "sequence_alignment_with_gaps_aa" after "query_sequence_aa"
cols = list(df.columns)
pos1 = cols.index("query_sequence_aa") + 1
cols.insert(pos1, cols.pop(cols.index("sequence_alignment_with_gaps_aa")))

# insert "sequence_alignment_aa_identical" after "aa_identical"
pos2 = cols.index("aa_identical") + 1
cols.insert(pos2, cols.pop(cols.index("sequence_alignment_aa_identical")))

df = df[cols]

# --- save with all logical values in uppercase ---
df = df.replace({True: "TRUE", False: "FALSE"})

# --- save the updated TSV ---
df.to_csv("add_aa_imgt.tsv", sep="\\t", index=False)

EOF

cp .command.log Translate_with_IMGT_gaps.log

    """
}
