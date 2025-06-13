#!/bin/bash

original="$1"     # original.fasta
aligned="$2"      # aligned.fasta
output="$3"       # restored.fasta

echo "Script has been entered"

# Créer une table de correspondance
grep '^>' "$original" | sed 's/^>//' > original_headers.txt
grep '^>' "$aligned" | sed 's/^>//' > aligned_headers.txt

paste aligned_headers.txt original_headers.txt > header_map.txt

echo "Temporary txt files have been made"

# Utiliser awk pour remplacer les headers dans le fichier aligné
awk '
  BEGIN {
    while ((getline < "header_map.txt") > 0) {
      split($0, a, "\t");
      map[a[1]] = a[2];
    }
  }
  /^>/ {
    header = substr($0, 2);
    if (header in map)
      print ">" map[header];
    else
      print $0;
    next;
  }
  { print }
' "$aligned" > "$output"

echo "Original aligned file :"
cat $aligned
echo "Restored headers file : "
cat $output