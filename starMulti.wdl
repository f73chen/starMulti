version 1.0

workflow starMulti {
input {
 Array[Pair[Pair[File, File], String]]+ inputFqsRgs
 String outputFileNamePrefix
}

scatter (fq_rg in inputFqsRgs) {
    File read1s    = fq_rg.left.left
    File read2s    = fq_rg.left.right
    String readgroups = fq_rg.right
}

call runStar { input: read1s = read1s, read2s = read2s, readgroups = readgroups, outputFileNamePrefix = outputFileNamePrefix }
call indexBam { input: inputBam = runStar.outputBam }

meta {
 author: "Peter Ruzanov & Jonathon Torchia"
 email: "peter.ruzanov@oicr.on.ca"
 description: "starMulti 1.0"
 dependencies: [
      {
        name: "star/2.7.3a",
        url: "https://github.com/alexdobin/STAR"
      },
      {
        name: "picard/2.19.2",
        url: "https://broadinstitute.github.io/picard/"
      }
    ]
}

output {
  File starBam          = runStar.outputBam
  File starChimeric     = runStar.outputChimeric
  File starIndex        = indexBam.outputBai
  File transcriptomeBam = runStar.transcriptomeBam
  File geneReadFile     = runStar.geneReads
 }
}

# ==========================================
#  TASK 1 of 2: run STAR aligner
# ==========================================
task runStar {
input {
  Array[File]+ read1s
  Array[File]+ read2s
  Array[String]+ readgroups
  String genome_index_dir = "$HG38_STAR_INDEX100_ROOT/"
  String outputFileNamePrefix
  String starSuffix = "Aligned.sortedByCoord.out"
  String transcriptomeSuffix = "Aligned.toTranscriptome.out"
  String chimericjunctionSuffix = "Chimeric.out"
  String genereadSuffix = "ReadsPerGene.out"
  String? addParam
  String modules = "star/2.7.3a hg38-star-index100/2.7.3a"
  Int uniqMAPQ = 255
  Int saSparsed = 2
  Int multiMax = -1
  Int chimSegmin = 12
  Int chimJunOvMin = 12
  Int chimOutJunFor = 1
  Int alignSJDBOvMin = 10
  Int alignMatGapMax = 100000
  Int alignIntMax = 100000
  Int chimMulmapScoRan = 3
  Int chimScoJunNonGTAG = -4
  Int chimMulmapNmax = 20
  Int chimNonchimScoDMin = 10
  Int peOvNbasesMin = 12
  Float peOvMMp = 0.1
  Int threads = 6
  Int jobMemory = 64
}

parameter_meta {
 read1s: "array of read1s"
 read2s: "array of read2s"
 readgroups: "array of readgroup lines"
 starSuffix: "Suffix for sorted file"
 transcriptomeSuffix: "Suffix for transcriptome-aligned file"
 chimericjunctionSuffix: "Suffix for chimeric junction file"
 genereadSuffix: "ReadsPerGene file suffix"
 addParam: "Additional STAR parameters"
 modules: "modules for running STAR"
 uniqMAPQ: "Score for unique mappers"
 saSparsed: "saSparsed parameter for STAR"
 multiMax: "multiMax parameter for STAR"
 chimSegmin: "minimum length of chimeric segment length"
 chimJunOvMin: "minimum overhang for a chimeric junction"
 chimOutJunFor: "formatting type for the Chimeric.out.junction  file"
 alignSJDBOvMin: "minimum overhang for annotated spliced alignments"
 alignMatGapMax: "maximum gap between two mates"
 alignIntMax: "maximum intron size"
 chimMulmapScoRan: "the score range for multi-mapping chimeras below the best chimeric score"
 chimScoJunNonGTAG: "penalty for a non-GTAG chimeric junction"
 chimMulmapNmax: "maximum number of chimeric multi-alignments"
 chimNonchimScoDMin: "to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value"
 peOvNbasesMin: "minimum number of overlap bases to trigger mates merging and realignment"
 peOvMMp: "maximum proportion of mismatched bases in the overlap area"
 threads: "Requested CPU threads"
 jobMemory: "Memory allocated for this job"
}

# missing --clip3pAdapterSeq $adaptors
command <<<
 STAR --twopassMode Basic \
      --genomeDir ~{genome_index_dir} \
      --readFilesIn ~{sep="," read1s} ~{sep="," read2s} \
      --readFilesCommand zcat \
      --outFilterIntronMotifs RemoveNoncanonical \
      --outFileNamePrefix ~{outputFileNamePrefix}. \
      --outSAMmultNmax ~{multiMax} \
      --outSAMattrRGline ~{sep=" , " readgroups} \
      --outSAMstrandField intronMotif \
      --outSAMmapqUnique  ~{uniqMAPQ} \
      --outSAMunmapped Within KeepPairs \
      --genomeSAsparseD ~{saSparsed} \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode TranscriptomeSAM GeneCounts \
      --chimSegmentMin ~{chimSegmin} \
      --chimJunctionOverhangMin ~{chimJunOvMin} \
      --chimOutJunctionFormat ~{chimOutJunFor} \
      --alignSJDBoverhangMin ~{alignSJDBOvMin} \
      --alignMatesGapMax ~{alignMatGapMax} \
      --alignIntronMax ~{alignIntMax} \
      --alignSJstitchMismatchNmax 5 -1 5 5 \
      --chimMultimapScoreRange ~{chimMulmapScoRan} \
      --chimScoreJunctionNonGTAG ~{chimScoJunNonGTAG} \
      --chimMultimapNmax ~{chimMulmapNmax} \
      --chimNonchimScoreDropMin ~{chimNonchimScoDMin} \
      --peOverlapNbasesMin ~{peOvNbasesMin} \
      --peOverlapMMp ~{peOvMMp} \
      --runThreadN ~{threads} ~{addParam}
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  cpu:     "~{threads}"
}


output {
 File outputBam        = "~{outputFileNamePrefix}.~{starSuffix}.bam"
 File outputChimeric   = "~{outputFileNamePrefix}.~{chimericjunctionSuffix}.junction"
 File transcriptomeBam = "~{outputFileNamePrefix}.~{transcriptomeSuffix}.bam"
 File geneReads        = "~{outputFileNamePrefix}.~{genereadSuffix}.tab"
}
}

# ==========================================
#  TASK 2 of 2: index bam file with picard
# ==========================================
task indexBam {
input {
	File   inputBam
  Int   jobMemory  = 12
  String? modules = "java/8 picard/2.19.2" 
}

parameter_meta {
 inputBam: "Input bam file"
 jobMemory: "Memory allocated indexing job"
 modules: "modules for running indexing job"
}

command <<<
 java -Xmx~{jobMemory-6}G -jar $PICARD_ROOT/picard.jar BuildBamIndex \
                              VALIDATION_STRINGENCY=LENIENT \
                              OUTPUT="~{basename(inputBam, '.bam')}.bai" \
                              INPUT=~{inputBam} 
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  File outputBai = "~{basename(inputBam, '.bam')}.bai"
}
}