#!/bin/zsh

# read mp1xxx.txt and generate image
input_file=$1
make run file=$input_file

# compare images
student_file="${input_file%.*}.png"
ref_file="mp1files/$student_file"

compare -fuzz 2% "$student_file" "$ref_file" ae.png
composite "$student_file" "$ref_file" -alpha off -compose difference rawdiff.png
convert rawdiff.png -level 0%,8% diff.png
convert +append "$ref_file" "$student_file" ae.png rawdiff.png diff.png "${student_file%.*}_result.png"

rm $student_file ae.png rawdiff.png diff.png
