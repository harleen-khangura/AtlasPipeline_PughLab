#!/bin/bash

#Enter email and API key associated with PEGR account and comma separated list of sample IDs and the respective TFs for analysis (e.g. 34544,34566 and GABPA,CTCF).
read -p "Enter User email:" USER_EMAIL
read -p "Enter PEGR API Key: " PEGR_API_KEY
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

    echo "Processing Sample ID: $SAMPLE_ID with transcription factor: $SAMPLE_TF"
    echo "$SAMPLE_ID" > "${SAMPLE_ID}.txt"

    # Run the Python script to generate the BAM file: EGC utility scripts required in folder for Python script to run 
    python3 generate_BAM_file_from_PEGR.py -f "${SAMPLE_ID}.txt" -p "$PEGR_API_KEY" -u "$USER_EMAIL" -b hg38

    # Determine the latest BAM file generated
    OUTPUT=$(ls *.bam | grep "$SAMPLE_ID")

    # Check if the BAM file exists
    if [ -f "$OUTPUT" ]; then
        NEW_BAM_NAME="${SAMPLE_ID}_${SAMPLE_TF}.bam"
        mv "$OUTPUT" "$NEW_BAM_NAME"
        echo "Renamed BAM file to: $NEW_BAM_NAME"

        samtools index "$NEW_BAM_NAME"
        echo "Indexed BAM file: $NEW_BAM_NAME"
    
        #If the BAM file exists, proceed with running BAM file through ChexMix with IgG master BAM control file. Must have control file (mergedoutput.bam) and Chexmix installed for this step
        java -Xmx10G -jar chexmix.v0.52.public.jar --geninfo hg38.info --expt "$NEW_BAM_NAME" --ctrl mergedoutput.bam --format BAM --out "${SAMPLE_ID}_${SAMPLE_TF}_chexmix" > "${SAMPLE_ID}_${SAMPLE_TF}_chexmix.out"


        # Run the Python script to expand ChexMix
        echo "Expanding ChexMix results..."

        EXPANDED_BED_FILE="/home/exouser/${SAMPLE_ID}${SAMPLE_TF}_chexmix_80bp.bed"
        java -jar ScriptManager-v0.15.jar coordinate-manipulation expand-bed -c=80 -o=$EXPANDED_BED_FILE "${SAMPLE_ID}_${SAMPLE_TF}_chexmix_experiment.bed"

       if [ -f "$EXPANDED_BED_FILE" ]; then
            echo "BED file created successfully: $EXPANDED_BED_FILE"

            #add two empty columns to bed file for formatting purposes
            FINAL_BED_FILE= "/home/exouser/final_${SAMPLE_ID}${SAMPLE_TF}_chexmix_80bp.bed"
            awk '{print $0 "\t.\t."}' $EXPANDED_BED_FILE  >  "/home/exouser/final_${SAMPLE_ID}${SAMPLE_TF}_chexmix_80bp.bed"       
            echo "Expanded Bed file rewritten"
            FINAL_BED_FILE = "/home/exouser/final_${SAMPLE_ID}${SAMPLE_TF}_chexmix_80bp.bed"

        # Run the Java command for sequence analysis
             echo "Running sequence analysis..."
             java -jar ScriptManager-v0.15.jar sequence-analysis fasta-extract --coord-header -o="${SAMPLE_ID}_chexmix.fasta" hg38.fa "/home/exouser/final_${SAMPLE_ID}${SAMPLE_TF}_chexmix_80bp.bed"
             echo "Processing for Sample ID ${SAMPLE_ID} completed."

        #move FASTA file to main home directory
         FASTA_NAME="${SAMPLE_ID}_chexmix.fasta"
         mv "$FASTA_NAME" ~/
         OUTPUT_DIR=~/"${SAMPLE_ID}${SAMPLE_TF}_memeresults"

        #Run MEME motif analysis through the CLI 
         echo "Running MEME analysis..."
         apptainer exec /home/exouser/meme.sif meme /home/exouser/"$FASTA_NAME"  -dna -oc "$OUTPUT_DIR" -nostatus -time 14400 -mod zoops -nmotifs 3 -minw 6 -maxw 50 -objfun classic -revcomp -markov_order 0
         echo "MEME analysis completed for Sample ID ${SAMPLE_ID}."

        #File organization so that all relevant files moved to corresponding Sample ID folder 
         mv "/home/exouser/final_${SAMPLE_ID}${SAMPLE_TF}_chexmix_80bp.bed" ~/"${SAMPLE_ID}${SAMPLE_TF}_memeresults" 
         mv "/home/exouser/${SAMPLE_ID}_${SAMPLE_TF}_chexmix_experiment.bed" ~/"${SAMPLE_ID}${SAMPLE_TF}_memeresults" 

        #run fimo 
         mv "$OUTPUT_DIR/meme.txt" ~/
         mv ~/meme.txt ~/"$SAMPLE_ID"meme.txt
        OUTPUT_DIR_2=~/"${SAMPLE_ID}${SAMPLE_TF}_motifvisualizations"
        apptainer exec /home/exouser/meme.sif fimo --oc "$OUTPUT_DIR_2" --verbosity 1 --bgfile --nrdb-- --thresh 1.0E-4 ~/"$SAMPLE_ID"meme.txt ~/"${SAMPLE_ID}_chexmix.fasta"


        else
            echo "Error: BED file not created for ${SAMPLE_ID}!"

        fi       

    else
        echo "BAM file for Sample ID $SAMPLE_ID not found!"
        echo "Available BAM files:"
        ls *.bam
    fi
done
