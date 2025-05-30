#!/bin/bash
for filename in inputs/$1/*$2*.toml; do
    if [[ $filename == *'3D'* ]]; then
        continue
    fi
    filename3D=${filename%%".toml"*}_3D${filename#*"${filename%%".toml"*}"}
    echo $filename3D
    cp $filename $filename3D
    echo -e "\nVolumetric = true" >> $filename3D
done
