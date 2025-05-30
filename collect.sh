#!/bin/bash
rm -rf outputs/$1/gamma
mkdir -p outputs/$1/gamma
for path in outputs/$1/*/$2/*.txt; do
    ln -f $path outputs/$1/gamma/$(basename ${path})
done
