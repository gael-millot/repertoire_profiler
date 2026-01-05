
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

def translate_with_gaps_and_mixed(nuc_seq):
    protein_seq = []
    mixed_positions = []  # store codon indices meeting new criteria

    seq_len = len(nuc_seq)
    valid_bases = set("ATGCatgc")

    # Track last mixed codon position
    last_mixed_pos = None

    for i in range(0, seq_len, 3):
        codon = nuc_seq[i:i+3]
        codon_index = i // 3 + 1  # codon position (1-based)

        # incomplete codon at the end
        if len(codon) < 3:
            protein_seq.append('-')
            continue

        # classify codon type
        is_all_dot = all(b == '.' for b in codon)
        is_all_dash = all(b == '-' for b in codon)
        has_valid = any(b in valid_bases for b in codon)
        has_invalid = any(b not in valid_bases for b in codon)
        is_mixed = has_valid and has_invalid
        no_valid = not has_valid  # codon fully without valid bases

        # --- translation logic ---
        if is_all_dot:
            protein_seq.append('.')
            continue
        elif is_all_dash:
            protein_seq.append('-')
            continue
        elif any(b not in valid_bases for b in codon):
            protein_seq.append('-')
        else:
            try:
                aa = str(Seq(codon.upper()).translate(to_stop=False))
            except Exception:
                aa = 'X'
            protein_seq.append(aa)

        # --- record logic ---
        if is_mixed:
            # condition 1: if we already saw a mixed codon before
            # condition 2: all codons between that and now contain no valid bases (dots or hyphens only)
            if last_mixed_pos is not None:
                mixed_positions.append(str(codon_index))
            last_mixed_pos = codon_index  # update last mixed position
        elif no_valid:
            # nothing to change, these can exist between mixed codons
            pass
        else:
            # if codon has valid bases only → break continuity
            last_mixed_pos = None

    protein = ''.join(protein_seq)
    mixed_str = ';'.join(mixed_positions)
    # return the aa seq and the mixed codon positions
    return protein, mixed_str 

# --- 1) import TSV file (header + single line) ---
df = pd.read_csv("${add_dotted_coord_ch}", sep="\\t")

# --- 2) get nucleotide sequence and 3) translate keeping dots and hyphens ---
def safe_translate(nuc_seq):
    # treat NA/None as missing: return missing outputs
    if pd.isna(nuc_seq):
        return (pd.NA, pd.NA)
    return translate_with_gaps_and_mixed(nuc_seq)

translated = df["sequence_alignment_with_gaps"].apply(safe_translate)
df["sequence_alignment_with_gaps_aa"] = translated.apply(lambda x: x[0]) # taken the first element
df["mixed_codon_positions"] = translated.apply(lambda x: x[1]) # taken the second element

# --- 5) compare AA sequences (ignore dots/hyphens) ---
def compare_ignore_gaps(seq1, seq2):
    if pd.isna(seq1) or pd.isna(seq2):
        return pd.NA

    s1 = seq1.replace('.', '').replace('-', '')
    s2 = seq2.replace('.', '').replace('-', '')
    return s1 == s2

df["sequence_alignment_aa_identical"] = df.apply(
    lambda x: compare_ignore_gaps(x["sequence_alignment_with_gaps_aa"], x["sequence_alignment_aa"]),
    axis=1
)

# --- 6 & 7) insert new columns in specific positions ---
# --- 6 & 7) insert new columns in specific positions ---
# move both "sequence_alignment_with_gaps_aa" and "mixed_codon_positions"
cols = list(df.columns)

# Find the insertion point (right after 'query_sequence_aa')
pos = cols.index("query_sequence_aa") + 1

# We remove and reinsert in reverse order so that the final order is preserved
for col in ["sequence_alignment_with_gaps_aa", "mixed_codon_positions"][::-1]:
    if col in cols:
        cols.insert(pos, cols.pop(cols.index(col)))

# insert "sequence_alignment_aa_identical" after "aa_identical"
pos2 = cols.index("aa_identical") + 1
cols.insert(pos2, cols.pop(cols.index("sequence_alignment_aa_identical")))

# Reorder the DataFrame
df = df[cols]

# --- save with all logical values in uppercase ---
df = df.replace({True: "TRUE", False: "FALSE"})

# --- save the updated TSV ---
df.to_csv("add_aa_imgt.tsv", sep="\\t", index=False)

EOF

cp .command.log Translate_with_IMGT_gaps.log

    """
}
