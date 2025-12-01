#!/bin/bash
start_time=$SECONDS
for filename in inputs/$1/*$2*.toml; do
    echo $filename
    go run . -i=$filename $3
done
elapsed_time=$((SECONDS - start_time))
minutes=$((elapsed_time / 60))
seconds=$((elapsed_time - minutes * 60))
echo "All done in $minutes m $seconds s"