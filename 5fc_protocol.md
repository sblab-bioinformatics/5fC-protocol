
# Basic pipeline to process sequencing runs from 5fC-enrichment method. 

## Aligning reads to the reference genome
### Data organization of raw `fastq` files.

Each replicate is stored in a different folder, named as `treatment_rep_X` or `control_rep_X`, where `X` in this case runs from 1 to 6. Each sequencing run was divide in lanes, such that the file names are of the form `SAMPLENAME_L???_R?_???.fastq.gz`, where `R?` is either `R1` or `R2` and indicate the paired ends.

Each job is sent to the compute cluster managed by a SLURM queuing system via `sbatch`.

## Adapater trimming

```bash
#! /usr/bin/env bash

declare -A new_names=( \
        ["treatment_rep_1"]="treatment_rep_1" \
        ["treatment_rep_2"]="treatment_rep_2" \
        ["treatment_rep_3"]="treatment_rep_3" \
        ["treatment_rep_4"]="treatment_rep_4" \
        ["treatment_rep_5"]="treatment_rep_5" \
        ["treatment_rep_6"]="treatment_rep_6" \
        ["control_rep_1"]="control_rep_1" \
        ["control_rep_2"]="control_rep_2" \
        ["control_rep_3"]="control_rep_3" \
        ["control_rep_4"]="control_rep_4" \
        ["control_rep_5"]="control_rep_5" \
        ["control_rep_6"]="control_rep_6" \
)

for folder in "${!new_names[@]}"
do
    pairs=()
    who=( ${new_names[$folder]}_L???_??_???.fastq.gz )
    for ff in ${folder}/${who[@]}
    do
        pairs+=("${ff%_??_???.*}")
    done
    sorted_unique_ids=($(echo "${pairs[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
    out_folder="trim_"${folder}

    if [ ! -d "${out_folder}" ]; then
        mkdir ${out_folder}
    fi

    for rr in  "${sorted_unique_ids[@]}"
    do
        cmd="trim_galore -q 10 --stringency 8 -o ${out_folder} --paired ${rr}_R1_001.fastq.gz ${rr}_R2_001.fastq.gz"
        echo ${cmd}
        sbatch -o ${rr}.log -J ${inp_folder}_${rr} --wrap "${cmd}" --mem 20G
    done
done
exit
```
### Alignment 

```bash
#! /usr/bin/env bash

declare -A new_names=( \
        ["treatment_rep_1"]="treatment_rep_1" \
        ["treatment_rep_2"]="treatment_rep_2" \
        ["treatment_rep_3"]="treatment_rep_3" \
        ["treatment_rep_4"]="treatment_rep_4" \
        ["treatment_rep_5"]="treatment_rep_5" \
        ["treatment_rep_6"]="treatment_rep_6" \
        ["control_rep_1"]="control_rep_1" \
        ["control_rep_2"]="control_rep_2" \
        ["control_rep_3"]="control_rep_3" \
        ["control_rep_4"]="control_rep_4" \
        ["control_rep_5"]="control_rep_5" \
        ["control_rep_6"]="control_rep_6" \
)

ref="~/mm9/Sequence/BWAIndex/genome.fa"

for folder in "${!new_names[@]}"
do
    inp_folder="trim_"${folder}
    bam_folder="bam_"${folder}
    who=( ${new_names[$folder]}_L???_??_???_val_?.fq.gz )
    pairs=()
    for ff in ${inp_folder}/${who[@]}
    do
        filename=$(basename "${ff}")
        short=${filename%_??_???_val_?.fq.gz}
        pairs+=("${short}")
    done
    sorted_unique_ids=($(echo "${pairs[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
    out_folder="bam_"${folder}

    if [ ! -d "${out_folder}" ]; then
        mkdir ${out_folder}
    fi

    for rr in  "${sorted_unique_ids[@]}"
    do
    cmd="bwa mem -M -t 8 $ref ${inp_folder}/${rr}_R1_001_val_1.fq.gz ${inp_folder}/${rr}_R2_001_val_2.fq.gz \
        | samtools sort -@8 - -o ${bam_folder}/${rr}.bam &&\
        samtools index ${bam_folder}/${rr}.bam"

    echo ${cmd}
    sbatch -o align_${rr}.log -J al_${inp_folder}_${rr} --wrap "${cmd}" --mem 30G
    done
done
```

### Merge files and filter out duplicates

```bash
#! /usr/bin/env bash

declare -A new_names=( \
        ["treatment_rep_1"]="treatment_rep_1" \
        ["treatment_rep_2"]="treatment_rep_2" \
        ["treatment_rep_3"]="treatment_rep_3" \
        ["treatment_rep_4"]="treatment_rep_4" \
        ["treatment_rep_5"]="treatment_rep_5" \
        ["treatment_rep_6"]="treatment_rep_6" \
        ["control_rep_1"]="control_rep_1" \
        ["control_rep_2"]="control_rep_2" \
        ["control_rep_3"]="control_rep_3" \
        ["control_rep_4"]="control_rep_4" \
        ["control_rep_5"]="control_rep_5" \
        ["control_rep_6"]="control_rep_6" \
)


picard_exe="~/sw/picard/picard-2.8.2/picard-2.8.2.jar"
whitelist="mm9-whitelist.bed"

for folder in "${!new_names[@]}"
do
    inp_folder="bam_"${folder}
    out_folder="bam_"${folder}
    who=( ${new_names[$folder]}_L???.bam )
    pairs=()
    for ff in ${inp_folder}/${who[@]}
    do
        pairs+=("${ff}")
    done

    name=${new_names[$folder]}
    cmd="samtools merge -f -@ 20 ${inp_folder}/merged_${name}.tmp.bam ${pairs[@]} && \
        java -Xmx3G -jar ${picard_exe} MarkDuplicates VALIDATION_STRINGENCY=SILENT I=${inp_folder}/merged_${name}.tmp.bam O=${inp_folder}/merged_${name}.bam M=${inp_folder}/${name}.markdup.txt &&\
        echo 'Filter bam files ' && \
        samtools view ${inp_folder}/merged_${name}.bam -f 3 -F 3840 -q 10 -b -L ${whitelist}  \
        | samtools sort -@ 20 -T /tmp/${name} -o ${inp_folder}/${name}.bam - && \
        echo ' Done sorting  ' && \
        samtools index ${inp_folder}/${name}.bam && \
        echo ' Done with all, cleaning up' && \
        rm -fr ${inp_folder}/merged_${name}.bam && \
        rm -fr ${inp_folder}/*_L00?.bam && \
        rm -fr ${inp_folder}/*_L00?.bam.bai && \
        echo 'Done' "

    echo ${cmd}
    sbatch -o merge_${name}.log -J m_${inp_folder}_${name} --wrap "${cmd}" --mem 60G
done
```

# Peak calling 

We use `MACS 2` for peak calling, pairing up the treatment and controls reads generated in the previous steps. 

```bash
#! /usr/bin/env bash

for rep in 1 2 3 4 5 6
do

out_folder=peak_calling_macs2_inp_${rep}

if [ ! -d "${out_folder}" ]; then
    mkdir ${out_folder}
fi
cd ${out_folder}

treatment_bam=bam_treatment_rep_"$rep"/treatment_rep_"$rep".bam
control_bam=bam_control_rep_"$rep"/control_rep_"$rep".bam

cmd="macs2 callpeak --keep-dup all --bdg -f BAMPE -t ../${treatment_bam} -c ../"$control_bam" -g mm -n peaks_rep"$rep"  "
echo "$cmd"
sbatch -o mcs2_${rep}.log -J mcs_${rep} --wrap "${cmd}" --mem 20G

cd ../
done

```

#### Extract peaks with 80% overlap accross replicates, 

Place all the `narrowPeak` files from `MACS2` in the same folder and then execute:

```bash
#!/usr/bin/env bash

function process() {
    # argument 1 --> associative array
    # argument 2 --> overlap value, must be [0,1], feed me well as I'm not testing it.

    eval "declare -A names="${1#*=}

    tableCat.py -i  ${names[@]} \
    | awk -v OFS="\t" '{print $1, $2, $3, $NF}' \
    | sort -k1,1 -k2,2n \
    | bedtools merge -c 4,4 -o distinct,count_distinct -i - \
    | gzip -f > 5fC_merged_peaks.bed.gz


    limit="$( echo ${#names[@]} | awk -v ov=$ov 'function ceil(x, y){y=int(x); return(x>y?y+1:y)} {print ceil($1*ov)}')"
    zcat 5fC_merged_peaks.bed.gz | awk -v limit="$limit" '($5>=limit) {print $0}' > 5fC_consensus_peaks.bed
    numb_consensus="$(wc -l 5fC_consensus_peaks.bed)"
    echo "Number of consensus sites using 80% overlap" "$numb_consensus"
    gzip -f 5fC_consensus_peaks.bed
}

declare -A names=(\
    ["rep_1"]="peaks_rep1_peaks.narrowPeak.gz" \
    ["rep_2"]="peaks_rep2_peaks.narrowPeak.gz" \
    ["rep_3"]="peaks_rep3_peaks.narrowPeak.gz" \
    ["rep_4"]="peaks_rep4_peaks.narrowPeak.gz" \
    ["rep_5"]="peaks_rep5_peaks.narrowPeak.gz" \
    ["rep_6"]="peaks_rep6_peaks.narrowPeak.gz" \
)

process "$(declare -p names)" 0.8
```

where tableCat.py can be found [here](scripts/tableCat.py).



