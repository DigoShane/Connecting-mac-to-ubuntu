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
            # Modify the specified lines using sed
            # sed -i 's/"lx" : 100.0,/"lx" : 10.0,/' "$file"
            # sed -i 's/"ly" : 100.0,/"ly" : 10.0,/' "$file"
            sed -i 's/"NN" : 10000/"NN" : 50000/' "$file"
            # sed -i 's/"read_in" : 0/"read_in" : 1/' "$file"
        else
            echo "File $file not found"
        fi
    else
        echo "Directory $dir not found"
    fi
done

