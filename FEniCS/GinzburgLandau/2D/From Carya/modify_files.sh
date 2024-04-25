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
            # Calculate the decimal expansion of i/100 up to 3 decimal places using bc
            decimal=$(echo "scale=3; $i/100" | bc)
            # Replace "H" : <some numbers> with "H" : $decimal using sed
            sed -i 's/"H" : [0-9]\+\.[0-9]\+/"H" : '"$decimal"'/' "$file" # sed is used for insertion, deletion, search and replace
            # Replace "read_in" : 0 with "read_in" : 1 using sed
            sed -i 's/"read_in" : 0/"read_in" : 1/' "$file"
            # Print the modified file content
            echo "Modified $file:"
            cat "$file"
            echo "============================================="
        else
            echo "File $file not found"
        fi
    else
        echo "Directory $dir not found"
    fi
done

