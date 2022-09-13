#! /usr/bin/env bash
set -euo pipefail

# By Thomas McCarthy @ BD

# Read in a datatable CSV such as the kind produced by BD Rhapsody Targeted Analysis Pipeline Version 1.10
# Convert to Matrix Market sparse matrix coordinate text file.
# This script takes filename ending ".csv" as its argument.

if test $# != 1 ; then
    echo "Number of arguments is not 1.
        Usage: bash ./convert_datatableCSV_to_MTX.sh <Input_CSV> 
          Input_CSV: a .csv file such as the kind produced by BD Rhapsody Targeted Analysis Pipeline Version 1.10.
        " 1>&2
    exit
fi

input="$1"
if test "${input: -4}" != ".csv" ; then
    echo "Input filename must end with '.csv'
        Usage: bash ./convert_datatableCSV_to_MTX.sh <Input_CSV> 
          Input_CSV: a .csv file such as the kind produced by BD Rhapsody Targeted Analysis Pipeline Version 1.10.
        " 1>&2
    exit
fi
# The output files go into a directory named after the input file.
out_directory=${input%.csv}
# If the output directory already exists, the script will error to avoid unintentional clobbering.
mkdir "$out_directory"



features_file="${out_directory}/features.tsv"
barcodes_file="${out_directory}/barcodes.tsv"
matrix_file="${out_directory}/matrix.mtx"
header_file="${out_directory}/temporary_header"
body_file="${out_directory}/temporary_body"
# Begin the header of the matrix file.
printf "%%%%MatrixMarket matrix coordinate integer general\n" > $header_file
printf "%%metadata_json: {}\n" >> $header_file

# Use one awk to strip out comment lines from the input.
# Use second awk to: add the features from the first non-comment line into the features file; add each barcode/cell-identifier to the barcodes file; store non-zero values in a temporary file; count the number of non-zero values in the matrix, which has to go at the top of the output.
# Antibody Capture features are differentiated from Gene Expression features by the string '|pAbO' at the end.
awk '$0 !~ /^#/ {print}' "$input" |
    awk -v FS=',' -v OFS=' ' -v features="$features_file" -v barcodes="$barcodes_file" -v headerfile="$header_file" -v bodyfile="$body_file" '
        BEGIN {count=0}
        NR == 1 {number_of_features=(NF - 1); 
            for(i=2; i<=NF; i++) {
                if ( $i ~ /\|pAbO$/ ) 
                    print $i "\t" $i "\t" "Antibody Capture" >>features 
                else 
                    print $i "\t" $i "\t" "Gene Expression" >>features
            }
            next
        }
        {
            print $1 "-1" >>barcodes; 
            for(i=2; i<NF; i++) {
                if ($i != 0) {
                    print (i - 1), (NR - 1), $i >>bodyfile;
                    count++
                };
            }
        }
        (NR % 1000) == 0 {print "Read " NR " lines."}
        END {number_of_barcodes=(NR - 1); print number_of_features, number_of_barcodes, count >>headerfile}
    '
# Concatenate the temporary header and body files into the matrix file.
cat "$header_file" "$body_file" > "$matrix_file"
rm "$header_file" "$body_file"
gzip "$features_file"
gzip "$barcodes_file"
gzip "$matrix_file"
