#!/bin/bash

nthreads = $1
genome_index = $2
chrom_sizes = $3

fastq1_rep1 = $4
fastq2_rep1 = $5
fastq1_rep2 = $6
fastq2_rep2 = $7

prefix_rep1 = $8
prefix_rep2 = $9
prefix = ${10}

outdir = ${11}
hic_analysis_path = ${12}

cd hic_analysis_path
wget https://raw.githubusercontent.com/4dn-dcic/docker-4dn-hic/refs/heads/master/scripts/run-merge-pairs.sh

cd outdir

# code adapted from https://github.com/4dn-dcic/docker-4dn-hic/tree/master

# Get pairs files for replicate 1

bwa mem -SP5M -t $nthreads $genome_index $fastq1_rep1 $fastq2_rep1 | samtools view -Shb - > $prefix_rep1.bam

samtools view -h $prefix_rep1.bam | {
	pairtools parse -c $chrom_sizes --add-columns mapq
} | {
	pairtools sort --nproc $nThreads --memory 32G --compress-program lz4c --tmpdir $outdir --output $prefix_rep1.sam.pairs.gz
}

pairtools dedup --mark-dups --output-dups - --output-unmapped - --output $prefix_rep1.marked.sam.pairs.gz $prefix_rep1.sam.pairs.gz

pairix $prefix_rep1.marked.sam.pairs.gz

pairtools stats --cmd-in 'pbgzip -dc -n '$nThreads'' -o $prefix_rep1.marked.pairs.stats $prefix_rep1.marked.sam.pairs.gz

pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' --output-rest $prefix_rep1.unmapped.sam.pairs.gz --output temp.gz $prefix_rep1.marked.sam.pairs.gz

pairtools split --output-pairs temp1.gz temp.gz

pairtools select 'True' --chrom-subset $chrom_sizes -o $prefix_rep1.dedup.pairs.gz temp1.gz

pairix $prefix_rep1.dedup.pairs.gz




# Get pairs files for replicate 2

bwa mem -SP5M -t $nthreads $genome_index $fastq1_rep2 $fastq2_rep2 | samtools view -Shb - > $prefix_rep2.bam

samtools view -h $prefix_rep2.bam | {
	pairtools parse -c $chrom_sizes --add-columns mapq
} | {
	pairtools sort --nproc $nThreads --memory 32G --compress-program lz4c --tmpdir $outdir --output $prefix_rep2.sam.pairs.gz
}

pairtools dedup --mark-dups --output-dups - --output-unmapped - --output $prefix_rep2.marked.sam.pairs.gz $prefix_rep2.sam.pairs.gz

pairix $prefix_rep2.marked.sam.pairs.gz

pairtools stats --cmd-in 'pbgzip -dc -n '$nThreads'' -o $prefix_rep2.marked.pairs.stats $prefix_rep2.marked.sam.pairs.gz

pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' --output-rest $prefix_rep2.unmapped.sam.pairs.gz --output temp.gz $prefix_rep2.marked.sam.pairs.gz

pairtools split --output-pairs temp1.gz temp.gz

pairtools select 'True' --chrom-subset $chrom_sizes -o $prefix_rep2.dedup.pairs.gz temp1.gz

pairix $prefix_rep2.dedup.pairs.gz


# Get fastq file from the two replicate pairs files

sh $hic_analysis_path/run_merge_pairs.sh $prefix.comb $prefix_rep1.dedup.pairs.gz $prefix_rep2.dedup.pairs.gz

cooler cload pairix -s 2 --assembly hg38 -p $nThreads $chrom_sizes:1000 $prefix.comb.pairs.gz $prefix.cool

cooler balance $prefix.cool

cooler zoomify -n $nThreads --balance --balance-args '--max-iters 500 --convergence-policy store_final' -r 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000 -o $prefix.mcool $prefix.cool







