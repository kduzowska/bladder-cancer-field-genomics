IN_DIR="/path/to/results/step4_pyclone_input"
OUT_DIR="/path/to/results/step5_pyclone_vi"

mkdir -p "$OUT_DIR"

for input in "$IN_DIR"/*_pyclone_input.tsv; do

    patient=$(basename "$input" "_pyclone_input.tsv")

    echo "Running $patient"

    pyclone-vi fit \
        -i "$input" \
        -o "$OUT_DIR/${patient}_pyclone.h5" \
        -c 40 \
        -d beta-binomial \
        -g 100 \
        -r 10 \
        -t 4 \
        --precision 200 \
        --seed 123 \
        > "$OUT_DIR/${patient}_fit.log" 2>&1

    pyclone-vi write-results-file \
        -i "$OUT_DIR/${patient}_pyclone.h5" \
        -o "$OUT_DIR/${patient}_pyclone_results.tsv"

done