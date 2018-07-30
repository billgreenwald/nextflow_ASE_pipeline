#!/usr/bin/env nextflow

//Configuration for the specific run//////////////////////////


outputDirectory="pipeline_out"
mode='paired' //should be either "single" or "paired" for single or paired end reads

vcfDir=file("testFiles/VCFs")
inputBam=file("testFiles/test.bam")
positionFile=file("testFiles/positions.tsv")


ramPerTry=4.GB
remapThreads=8

/////////////////////////////////////////////////////////////

//Constants for the pipeline

//scripts
WaspDir=file("bin/WASP/")
getFn=file("scripts/getFileName.py")
BWA=file("bin/bwa-0.7.15/bwa")
samtools=file("bin/samtools-1.4/samtools")
processPileup=file("scripts/Process_Pileup.py")


//files
referenceGenome=file("files/hg19_sorted.fa")

process makeSNPDir{
	executor 'sge'
	cpus 1
	memory { 2 * ramPerTry}
    queue "short"
    clusterOptions "-l short -S /bin/bash -l h_vmem=${task.memory.toGiga()}G"
    penv 'smp'
    
	output:
	file "SNPs/" into SNPdir

	script:
    """
    mkdir SNPs
    for i in `seq 1 22`; do
        awk -v OFS="\\t" '{if (substr(\$1,1,1)!="#") print \$2,\$4,\$5}' $vcfDir/chr\$i.vcf > SNPs/chr\$i.snps.txt
        gzip SNPs/chr\$i.snps.txt
    done
    """
}

process findIntersectingSnps{
	executor 'sge'
	cpus 1
	memory { 2 * ramPerTry}
    queue "short"
    clusterOptions "-l short -S /bin/bash -l h_vmem=${task.memory.toGiga()}G"
    penv 'smp'
    
    input:
    file SNPdir
        
	output:
    file "remap.fq1.gz" into toRemap1
    file "remap.fq2.gz" into toRemap2
    file "keep.bam" into keepBam 
    file "to.remap.bam" into toRemapBam

	script:
    if (mode == "paired")
    """
    python $WaspDir/mapping/find_intersecting_snps.py --is_paired_end --is_sorted $inputBam --snp_dir $SNPdir --output_dir ./
    fileName=\$(python $getFn $inputBam)
    mv \$fileName.keep.bam keep.bam
    mv \$fileName.remap.fq1.gz remap.fq1.gz
    mv \$fileName.remap.fq2.gz remap.fq2.gz
    mv \$fileName.to.remap.bam to.remap.bam
    """
    else if (mode == "single")
    """
    python $WaspDir/mapping/find_intersecting_snps.py --is_sorted $inputBam --snp_dir $SNPdir --output_dir ./
    fileName=\$(python $getFn $inputBam)
    mv \$fileName.keep.bam keep.bam
    mv \$fileName.remap.fq1.gz remap.fq1.gz
    mv \$fileName.toRemap.bam toRemap.bam
    echo "File created to make nextflow pipeline run" > $remap.fq2.gz
    """    
}

process remap{
	executor 'sge'
	cpus remapThreads
	memory { 2 * ramPerTry}
    queue "short"
    clusterOptions "-l short -S /bin/bash -l h_vmem=${task.memory.toGiga()}G"
    penv 'smp'
    
    input:
    file toRemap1
    file toRemap2
        
	output:
	file "Aligned.out.bam" into remapped

	script:
    if (mode=="paired")
    """
    $BWA mem -t $remapThreads $referenceGenome $toRemap1 $toRemap2 > Aligned.out.sam

    $samtools view -b Aligned.out.sam > Aligned.out.bam
    """ 
    
    else if (mode=="single")  
    """
    $BWA mem -t $remapThreads $referenceGenome $toRemap1 > $Aligned.out.sam

    $samtools view -b Aligned.out.sam > Aligned.out.bam
    """ 
}

process sortRmNoSup{
	executor 'sge'
	cpus 1
	memory { 2 * ramPerTry}
    queue "short"
    clusterOptions "-l short -S /bin/bash -l h_vmem=${task.memory.toGiga()}G"
    penv 'smp'
    
    input:
    file remapped
    file keepBam
    file toRemapBam
        
	output:
	file "remapped.merged.NoSup.sorted.bam" into remapMerge
    file "remapped.merged.NoSup.sorted.bam.bai" into remapMergeBai

	script:
    """
    python $WaspDir/mapping/filter_remapped_reads.py $toRemapBam  $remapped remappedOutput.bam 

    $samtools merge -f remapped.merged.bam remappedOutput.bam $keepBam
    $samtools view -bF 2048 remapped.merged.bam > remapped.merged.NoSup.bam
    $samtools sort  remapped.merged.NoSup.bam > remapped.merged.NoSup.sorted.bam
    $samtools index remapped.merged.NoSup.sorted.bam
    """ 
}

process rmDupe{
    publishDir outputDirectory, mode: "copy"
	executor 'sge'
	cpus 1
	memory { 2 * ramPerTry}
    queue "short"
    clusterOptions "-l short -S /bin/bash -l h_vmem=${task.memory.toGiga()}G"
    penv 'smp'
    
    input:
    file remapMerge
    file remapMergeBai
        
	output:
	file "remapped.merged.NoSup.sorted.rm_dup.sorted.bam" into waspOut

	script:
    if (mode=="paired")
    """
    python $WaspDir/mapping/rmdup_pe.py $remapMerge remapped.merged.NoSup.sorted.rm_dup.bam
    $samtools sort remapped.merged.NoSup.sorted.rm_dup.bam > remapped.merged.NoSup.sorted.rm_dup.sorted.bam
    """
    
    else if (mode == "single")
    """
    python $WaspDir/mapping/rmdup.py $remapMerge remapped.merged.NoSup.sorted.rm_dup.bam
    $samtools sort remapped.merged.NoSup.sorted.rm_dup.bam > remapped.merged.NoSup.sorted.rm_dup.sorted.bam
    """
}

process pileup{
    publishDir outputDirectory, mode: "copy"
	executor 'sge'
	cpus 1
	memory { 2 * ramPerTry}
    queue "short"
    clusterOptions "-l short -S /bin/bash -l h_vmem=${task.memory.toGiga()}G"
    penv 'smp'

    
    input:
    file waspOut
        
	output:
	file "pileup.tsv" into pileup
    file "processedPileup.tsv" into pileup2
    
	script:
    """
    $samtools mpileup -l $positionFile -f $referenceGenome $waspOut > pileup.tsv
    python $processPileup pileup.tsv processedPileup.tsv
    """
}