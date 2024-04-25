#!/bin/bash

# Iterate over directories H0 to H100
for ((i=0; i<=100; i++)); do
    dir="H$i"
    if [ -d "$dir" ]; then
        # Check if the directory exists
        file="$dir/Input$i.json"
        if [ -f "$file" ]; then
            # Check if the file exists
            echo "Modifying $file"
            cat "$file"
            echo "============================================="
        else
            echo "File $file not found"
        fi
    else
        echo "Directory $dir not found"
    fi
done

