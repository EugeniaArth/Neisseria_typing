#!/bin/bash

# Update paths for  GATK
export PATH="gatk-4.3.0.0/:$PATH"
export PATH="picard.jar:$PATH"

# Define reference files and directories
REFERENCE="/reference/FA19.fasta"
REF_ANNOTATION="/reference/FA19.gff3"
ANNOVAR_DB="/GATK4/annovar_db/"
RESULTS_FILE="/GATK4/calling_results.txt"


# Base directories
MAPPED_DIR="/GATK4/mapped/"
CALLING_DIR="/GATK4/calling/"
LOG_DIR="/GATK4/calling/logs/"

# Create log directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Initialize results file with header
echo -e "Sample Name\tCountry\tCalling" > "$RESULTS_FILE"

# Create sequence dictionary and index for both references if dont have
samtools faidx $REFERENCE
java -jar picard.jar CreateSequenceDictionary R=$REFERENCE O=${REFERENCE%.fasta}.dict



# Process each country directory in the mapped directory
for COUNTRY_DIR in $MAPPED_DIR*/; do
    COUNTRY=$(basename "$COUNTRY_DIR")
    echo "Processing country: $COUNTRY"

    # Create necessary directories for this country
    BAM_DIR=${CALLING_DIR}/${COUNTRY}/bam
    VCF_DIR=${CALLING_DIR}/${COUNTRY}/vcf
    ANNOTATED_VCF_DIR=${CALLING_DIR}/${COUNTRY}/vcf_annotated
    mkdir -p $BAM_DIR $VCF_DIR $ANNOTATED_VCF_DIR

    # Process each sorted BAM file in the country's directory
    for BAM_FILE in ${COUNTRY_DIR}*.sorted.bam; do
        echo "Processing file: $BAM_FILE"

        # Extract base name from BAM file
        BASE_NAME=$(basename "$BAM_FILE" .sorted.bam)
        LOG_FILE="$LOG_DIR/${BASE_NAME}_calling.log"


        # Add read groups, sort, and index
        java -jar picard.jar AddOrReplaceReadGroups \
            -I $BAM_FILE \
            -O $BAM_DIR/${BASE_NAME}.rgs.bam \
            -RGID $BASE_NAME -RGLB miseq -RGSM $BASE_NAME -RGPL illumina -RGPU unit1
        samtools index $BAM_DIR/${BASE_NAME}.rgs.bam

        # Sort and index the BAM file with read groups
        samtools sort $BAM_DIR/${BASE_NAME}.rgs.bam -o $BAM_DIR/${BASE_NAME}.rgs.sorted.bam
        samtools index $BAM_DIR/${BASE_NAME}.rgs.sorted.bam

        # Mark duplicates and index
        java -jar picard.jar MarkDuplicates \
            -I $BAM_DIR/${BASE_NAME}.rgs.sorted.bam \
            -O $BAM_DIR/${BASE_NAME}.marked_dup.bam \
            -M $BAM_DIR/${BASE_NAME}.dup_metrics.txt
        samtools index $BAM_DIR/${BASE_NAME}.marked_dup.bam

        # Call variants using GATK HaplotypeCaller
        gatk HaplotypeCaller --native-pair-hmm-threads 8 \
            -ploidy 1 \
            -R $REFERENCE \
            -I $BAM_DIR/${BASE_NAME}.marked_dup.bam  \
            -O $VCF_DIR/${BASE_NAME}.vcf

            # Extract SNPs & Indels
          gatk SelectVariants \
                    -R $REFERENCE \
                    -V $VCF_DIR/${BASE_NAME}.vcf \
                    --select-type SNP \
                    -O $VCF_DIR/${BASE_NAME}.snp.vcf
          echo SNP extracted for ${BASE_NAME}

          gatk SelectVariants \
                    -R $REFERENCE \
                    -V $VCF_DIR/${BASE_NAME}.vcf \
                    --select-type INDEL \
                    -O $VCF_DIR/${BASE_NAME}.indel.vcf
          echo INDEL extracted for ${BASE_NAME}

        # Filter variants based on quality metrics
        # SNP
        gatk --java-options "-Xmx4G" VariantFiltration -R $REFERENCE -V $VCF_DIR/${BASE_NAME}.snp.vcf \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O $VCF_DIR/${BASE_NAME}.snp.filtered.vcf

        echo SNP filtration is performed for ${BASE_NAME}

        # indels
        gatk --java-options "-Xmx4G" VariantFiltration -R $REFERENCE -V $VCF_DIR/${BASE_NAME}.indel.vcf \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O $VCF_DIR/${BASE_NAME}.indel.filtered.vcf

        echo Indels filtration is performed for ${BASE_NAME}

     # убираем отфильтрованные / Exclude Filtered Variants
       gatk --java-options "-Xmx4G" SelectVariants --exclude-filtered -V $VCF_DIR/${BASE_NAME}.snp.filtered.vcf   \
       -O $VCF_DIR/${BASE_NAME}.snp.passed.vcf

       echo variants filtered for ${BASE_NAME}

       gatk --java-options "-Xmx4G" SelectVariants --exclude-filtered -V $VCF_DIR/${BASE_NAME}.indel.filtered.vcf  \
        -O $VCF_DIR/${BASE_NAME}.indel.passed.vcf
       echo variants filtered for ${BASE_NAME}



    # name genes by simply intersecting with reference
    bedtools intersect -wb -a $VCF_DIR/${BASE_NAME}.snp.passed.vcf  -b $REF_ANNOTATION  -header >  $VCF_DIR/${BASE_NAME}.snp.named.vcf
    bedtools intersect -wb -a $VCF_DIR/${BASE_NAME}.indel.passed.vcf -b $REF_ANNOTATION -header >  $VCF_DIR/${BASE_NAME}.indel.named.vcf
    echo bedtools intersect is finished for ${BASE_NAME}


    # Next step is annotation of variants using Annovar

    ## Need to change directory as annovar will create output in the working directory
    cd $VCF_DIR



    #Annotate SNPs and Predict Effects
    perl  /annovar/table_annovar.pl  $VCF_DIR/${BASE_NAME}.snp.named.vcf $ANNOVAR_DB --outfile $ANNOTATED_VCF_DIR/${BASE_NAME}.snp --protocol \
    refGene --operation g --buildver  FA19 --vcfinput
    echo SNP annotation is performed for ${BASE_NAME}

    perl /annovar/table_annovar.pl  $VCF_DIR/${BASE_NAME}.indel.named.vcf $ANNOVAR_DB --outfile $ANNOTATED_VCF_DIR/${BASE_NAME}.indel --protocol \
    refGene --operation g --buildver  FA19 --vcfinput

    echo INDELS annotation is performed for ${BASE_NAME}


      done
  done
