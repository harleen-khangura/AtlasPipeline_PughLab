
#index files
directory="/media/volume/atlas-data"

for bam_file in *.bam; do

    bai_file="${bam_file}.bai"

    if [ -f "$bai_file" ]; then
    echo "Index for $bam_file already exists, skipping indexing."
    continue
    fi

    echo "Indexing $bam_file..."
    samtools index "$bam_file"

    if [ $? -eq 0 ]; then
        echo "Successfully indexed $bam_file"
    else
        echo "Failed to index $bam_file" >&2
    fi
done

#use only bai files
input_dir="/media/volume/atlas-data/"
bai_files=$(find "$input_dir" -type f -name "*.bai")

if [ -z "$bai_files" ]; then
    echo "No .bai files found in the directory."
else
    echo "List of .bai files found:"
    echo "$bai_files"
fi


# Run BAM correlation 
#did not update outpute file (kind of fixed)
files=("38188_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam" "34031_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam"  
"34055_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam"
"37449_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam"
"33925_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam"
"38471_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam"
"38631_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam")
                                                
num_files=${#files[@]}                          
                          
for ((i=0; i<num_files; i++)); do
    for ((j=i+1; j<num_files; j++)); do     
        file1=${files[$i]} 
        file2=${files[$j]}
        echo "Pair: $file1 and $file2"
        java -jar ScriptManager-v0.15.jar bam-statistics bam-correlation $file1 $file2  --output=correlation_results_$file1_and_$file2.txt
    done                              
done 

#put values into matrix 

matrix_file="final_correlation_matrix.txt"

row1="entryid 33925	34031	34055	37449	38188	38471	38631"
row2="33925	—	0.0712017831175956	0.0781389878383116	0.143056986170486	0.00401975197164926	0.101386599711553	0.0789347626676286"
row3="34031	0.0712017831175956	—	0.0997001808327295	0.0886377703761014	0.00399369995843822	0.0672621432423906	0.0502736566032051"
row4="34055	0.0781389878383116	0.0997001808327295	—	0.0914557991457	0.00387518157954927	0.0705976980945756	0.0509729380861669"
row5="37449	0.143056986170486	0.0886377703761014	0.0914557991457	—	0.00773068128505752	0.249047356982112	0.202890900410771"
row6="38188	0.00401975197164926	0.00399369995843822	0.00387518157954927	0.00773068128505752	—	0.00581024959688649	0.0034150760777767"
row7="38471	0.101386599711553	0.0672621432423906	0.0705976980945756	0.249047356982112	0.00581024959688649	—	0.131136967805419"
row8="38631	0.0789347626676286	0.0502736566032051	0.0509729380861669	0.202890900410771	0.0034150760777767	0.131136967805419	—"

echo -n "" > "$matrix_file" 

echo "$row1" >> "$matrix_file"
echo "$row2" >> "$matrix_file"
echo "$row3" >> "$matrix_file"
echo "$row4" >> "$matrix_file"
echo "$row5" >> "$matrix_file"
echo "$row6" >> "$matrix_file"
echo "$row7" >> "$matrix_file"
echo "$row8" >> "$matrix_file"

echo "Matrix file created and populated: $matrix_file"


#merge files excluding outlier 
samtools merge -f mergedoutput.bam 33925_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam 34031_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam 34055_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam 37449_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam 38471_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam 38631_IgG_i5006_K562_-_IMDM_-_BX_hg38.bam

