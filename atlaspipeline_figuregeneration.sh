#!/bin/bash


#read -p "Enter User email:" USER_EMAIL
#read -p "Enter PEGR API Key: " PEGR_API_KEY
read -p "Enter a comma-separated list of sample IDs: " SAMPLE_IDS
read -p "Enter a comma-separated list of transcription factors for each sample: " SAMPLE_TFS

echo "Sample IDs and transcription factors entered for processing: $SAMPLE_IDS and $SAMPLE_TFS" 

# Convert SAMPLE_IDS and SAMPLE_TFS into arrays
IFS=',' read -r -a SAMPLE_ID_ARRAY <<< "$SAMPLE_IDS"
IFS=',' read -r -a SAMPLE_TF_ARRAY <<< "$SAMPLE_TFS"

# Check if both arrays are the same length
if [ "${#SAMPLE_ID_ARRAY[@]}" -ne "${#SAMPLE_TF_ARRAY[@]}" ]; then
    echo "Error: The number of sample IDs does not match the number of transcription factors."
    exit 1
fi

# Process each sample ID with its corresponding transcription factor
for i in "${!SAMPLE_ID_ARRAY[@]}"; do
    # Trim whitespace
    SAMPLE_ID=$(echo "${SAMPLE_ID_ARRAY[$i]}" | xargs)
    SAMPLE_TF=$(echo "${SAMPLE_TF_ARRAY[$i]}" | xargs)

    OUTPUT_DIR_2=~/"${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations"

    #Separate FIMO output GFF into three separate files containing lines for individual motifs: writefimomotifs script can be found on Github 
    python3 ~/writefimomotifs.py "$OUTPUT_DIR_2"/fimo.gff "$OUTPUT_DIR_2"
    mv ~/vol/"${SAMPLE_ID}_${SAMPLE_TF}.bam" ~/vol/"${SAMPLE_ID}_${SAMPLE_TF}.bam.bai"  "$OUTPUT_DIR_2"
    FASTA_NAME= /home/exouser/"${SAMPLE_ID}_chexmix.fasta"

        #Generate four colour plot (convert gff to bed and fasta w/ Chexmix peaks as reference genome for each of the three motifs and move to Motif Visualizations folder)
for i in {1..3}
        do
            input_file="${OUTPUT_DIR_2}/meme_${i}_fimo.gff" 
            output_file="${OUTPUT_DIR_2}/meme_${i}_fimo.bed" 
            output_fasta_file="${OUTPUT_DIR_2}/meme_${i}_fimo.fasta" 
            java -jar ScriptManager.jar coordinate-manipulation gff-to-bed -o="${output_file}" "${input_file}"
           echo " java -jar ScriptManager-v0.15.jar sequence-analysis fasta-extract --coord-header -o /home/exouser/"${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/meme_${i}_fimo.fasta" ~/"${SAMPLE_ID}_chexmix.fasta" "${output_file}""
           java -jar ScriptManager-v0.15.jar sequence-analysis fasta-extract --coord-header -o /home/exouser/"${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/meme_${i}_fimo.fasta" ~/"${SAMPLE_ID}_chexmix.fasta"  "${output_file}"  
            java -jar ScriptManager-v0.15.jar figure-generation four-color -o="${OUTPUT_DIR_2}/${SAMPLE_ID}_motif${i}.png" -x=1 -y=1 "${output_fasta_file}"
        done
        
        #Composite plot generation 
            # Make new fimo bed with newbedforcomposite Python script to adjust chromosome start and end positions
            #Expand bed by 1000 bp
            #Tagpileup + figure generation 
for i in {1..3}
        do
           OUTPUT_DIR_2=~/"${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations"
            python3 ~/newbedforcomposite.py "${OUTPUT_DIR_2}/meme_${i}_fimo.bed" "${OUTPUT_DIR_2}/meme${i}fimo_tagpileup.bed" 
            java -jar ScriptManager.jar coordinate-manipulation expand-bed -c=1000 -o="/home/exouser/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/meme${i}fimo_tagpileup_1000bp.bed"  "/home/exouser/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/meme${i}fimo_tagpileup.bed" 
            java -jar ScriptManager.jar read-analysis tag-pileup     -a   -M="${OUTPUT_DIR_2}/meme${i}_matrix"  -o="${OUTPUT_DIR_2}/meme${i}_composite"  -b=1 ~/"${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/meme${i}fimo_tagpileup_1000bp.bed" ~/"${SAMPLE_ID}${SAMPLE_TF}_motifvisualization/"${SAMPLE_ID}_${SAMPLE_TF}.bam
            java -jar ScriptManager.jar figure-generation composite-plot -o="${OUTPUT_DIR_2}/meme${i}plot.png" -l "${OUTPUT_DIR_2}/meme${i}_composite" 

  
        done

#Heat map generation 
#Sort bed by window occupancy and tag pileup 
# Generate forward and reverse strand heat map and merge 
#Gzip heat map or can add commented line to only produce gzipped heat map instead of both 

for i in {1..3}
        do
           awk '{print $1"\t"$2"\t"$3"\tmeme1_row(" NR ")\t"$5"\t"$6}' "/home/exouser/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/meme${i}fimo_tagpileup_1000bp.bed" > "/home/exouser/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/labelled_meme${i}fimo_tagpileup_1000bp.bed"
            java -jar ScriptManager.jar read-analysis tag-pileup     -a  --combined  -M="/home/exouser/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/heatmapmeme${i}_matrix"  -o="/home/exouser/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/heatmapmeme${i}_composite"  -b=1 "/home/exouser/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/labelled_meme${i}fimo_tagpileup_1000bp.bed" "~/"${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/"${SAMPLE_ID}_${SAMPLE_TF}.bam"
            java -jar ScriptManager.jar coordinate-manipulation sort-bed -o "~/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/sortedbed_forheatmap" "~/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/labelled_meme${i}fimo_tagpileup_1000bp.bed" "heatmapmeme${i}_matrix_combined.cdt"
            java -jar ScriptManager.jar read-analysis tag-pileup     -a   -M "~/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/finalheatmapmeme${i}_matrix"   -o "~/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/finalheatmapmeme${i}_composite"  -b=1 "~/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/sortedbed_forheatmap.bed" "~/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/${SAMPLE_ID}_${SAMPLE_TF}.bam"
            java -jar ScriptManager.jar figure-generation heatmap --blue -o="/home/exouser/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/heatmap_sense_meme${i}.png" -p=0.95 "~/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/finalheatmapmeme${i}_matrix_sense.cdt"
            java -jar ScriptManager.jar figure-generation heatmap --red -o="/home/exouser/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/heatmap_anti_meme${i}.png" -p=0.95 "~/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/finalheatmapmeme${i}_matrix_anti.cdt"
            java -jar ScriptManager.jar figure-generation merge-heatmap "/home/exouser/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/heatmap_sense_meme${i}.png" "/home/exouser/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/heatmap_anti_meme${i}.png" -o="/home/exouser/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/heatmap_combined_meme${i}.png" 
            #gzip -c ~/39382ZNF276_motifvisualizations/heatmaptest.png > ~/39382ZNF276_motifvisualizations/heatmaptest.png.gz
            gzip "/home/exouser/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/heatmap_combined_meme${i}.png" > "/home/exouser/${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations/heatmap_combined_meme${i}.png.gz"

        done
    
#Move relevant files to Sample's Motif Visualizations folder 
mv  ~/"${SAMPLE_ID}_chexmix.fasta.fai"  ~/"${SAMPLE_ID}_chexmix.fasta" ~/"${SAMPLE_ID}meme.txt" ~/"${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations"    
done 
