#!/bin/bash
for filename in inputs/$1/*$2*.toml; do
    echo $filename
done
