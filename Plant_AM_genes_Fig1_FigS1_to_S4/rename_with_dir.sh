#!/bin/bash
shopt -s globstar  ##globstar will let us match files recursively
files=( /home/pers/juliana.almario/gemoma/**/predicted_proteins.fasta )  ##Array containing matched files, mention where to search and what files here
for i in "${files[@]}"; do 
    d="${i%/*}"  ##Parameter expansion, gets the path upto the parent directory
    d_="${d##*/}"  ##gets the name of parent directory
    f="${i##*/}"  ##gets the file name
     echo mv "$i" "$d"/"${d_}""_""$f"  ##renaming, remove echo after confirming what will be changed and you are good
