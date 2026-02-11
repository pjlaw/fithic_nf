


def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
      nextflow run fithic.nf --hicPath /path/to/hic/files --resolution BP --outDir /path/to/output

    Mandatory arguments:
      --hicPath /path/to/hic/files      Folder with all the hic files. All *hic files will be analysed
      --resolution BP                   Resolution to run the analysis at
      --outDir /path/to/output          Output folder
	
	
    Optional arguments:
      --fdr Q           Pre-merge Q-value filter (only merge those interactions with FDR<Q). Default: 0.01
	  --washUFilter Q   Filter interactions for the washU bedfile. Default:1e-10
	  --distFilter Q    Remove interactions with distance < Q. Default:10000

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}


/*
 * load HiC files
 */ 

Channel
    .fromPath("${params.hicPath}/*.hic")
    .set {hicFiles_ch}
    
//half the resolution for half bin size
//half_bin=params.resolution.intdiv(2)

Channel.of(1..22, 'X')
	.set {chromosomes_ch}


//convert hic to contactMap
process dumpHic {
    tag "${hic_file}_dump"
    publishDir "${params.outDir}/dump/", mode: "copy", overwrite: true
    
    
	input:
		file hic_file from hicFiles_ch
		each chrom from chromosomes_ch
	
	//filename="${hic_file.baseName}_chr${chrom}_${params.resolution}"
	
	output:		
		tuple val(chrom), env(filename), env(experiment_name), file("*.contactCounts.gz") into contacts_fithic_ch
		tuple val(chrom), env(filename), env(experiment_name), env(normalisation), file("*.norm") into normalisation_ch
	    tuple val(chrom), env(filename), env(experiment_name), file("${params.build}_chr${chrom}_${params.resolution}_fithic_fragmentsfile.gz") into (mapbias_ch, map_fithic_ch)
	    
	//technically already loaded to run nextflow
	module 'java'
	
	shell:
	'''
	#create mappibility files
	if [ -f !{params.mappabilityPath}/!{params.build}_chr!{chrom}_!{params.resolution}_fithic_fragmentsfile.gz ]; then
		ln -s !{params.mappabilityPath}/!{params.build}_chr!{chrom}_!{params.resolution}_fithic_fragmentsfile.gz
	else
		#NOTE: modified version that one processes one chromosome
		createFitHiCFragments-fixedsize_modified.py --chrLens !{params.chrLensPath} --resolution !{params.resolution} --chrom chr!{chrom} --outFile !{params.build}_chr!{chrom}_!{params.resolution}_fithic_fragmentsfile.gz
	fi
	
	experiment_name=$(basename !{hic_file} .hic)
	chromosome="chr!{chrom}"
	filename="${experiment_name}_${chromosome}_!{params.resolution}"

	export _JAVA_OPTIONS="-Xmx8g -Xms8g"  
	#dump out the interactions
	java -jar !{params.juicerTools} dump observed NONE !{hic_file} $chromosome $chromosome BP !{params.resolution} ${filename}.contacts.raw
	#shift to midpoints
	#awk -v half=${half_bin} '{print $1"\t"$2-half"\t"$3"\t"$4-half"\t"$5}' "${filename}.contacts" | gzip > "${filename}_midshifted.contacts"
	
	#from createFitHiCContacts-hic.sh 
	awk -v chr1=$chromosome -v chr2=$chromosome -v r=!{params.resolution} '{OFS="\t"} { print chr1,$1-r/2,chr2,$2-r/2,$3}' ${filename}.contacts.raw | gzip > ${filename}.contactCounts.gz 
	
	#get the normalisations from juicer, rather than the fithic script
    function try()
    {
        [[ $- = *e* ]]; SAVED_OPT_E=$?
        set +e
    }
    
    function catch()
    {
        export exception_code=$?
        (( $SAVED_OPT_E )) && set +e
        return $exception_code
    }
    
    normalisation="KR"
    try; (
     java -jar !{params.juicerTools} dump norm ${normalisation} !{hic_file} $chromosome BP !{params.resolution} ${filename}_${normalisation}.norm
    )
    catch || {
        case $exception_code in
        9)
            rm ${filename}_${normalisation}.norm
            normalisation="VC"
            java -jar !{params.juicerTools} dump norm ${normalisation} !{hic_file} $chromosome BP !{params.resolution} ${filename}_${normalisation}.norm
            ;;
        *)
            echo "unknown error"
            ;;
        esac
    }
	
	
	'''
}

// using the juicer norm dump rather than HiCKRy.py
process computeBiases {
    tag "${filename}_bias"
    publishDir "${params.outDir}/bias/", mode: "copy", overwrite: true

    input:
		tuple val(chrom), val(filename), val(experiment_name), file('mappabilityFile.gz'),  val(normalisation), file('normfile.txt') from mapbias_ch.join(normalisation_ch, by: [0,1,2])
		
	output:
		tuple val(chrom), val(filename), val(experiment_name), file("${filename}_${normalisation}_bias.gz") into bias_fithic_ch
	
	script:
	"""
	#replace NaN with -1
	sed -i 's/NaN/-1/g' normfile.txt
	#extract chromsome and midpoint
	zcat mappabilityFile.gz | cut -f1,3 > tempfile
	#add the normalised values from juicer
	paste tempfile normfile.txt | gzip > ${filename}_${normalisation}_bias.gz
	"""
}

 //val(chrom), val(filename), file("${filename}_fithic.contactCounts.gz"),file("${filename}_bias.gz")
contacts_fithic_ch.join(bias_fithic_ch, by:[0,1,2])
	.set{ contacts_bias_ch }

process runFitHic {
    tag "${filename}_fithic"
    publishDir "${params.outDir}/fithic/", mode: "copy", overwrite: true
    errorStrategy 'retry'
	maxRetries = 2
    cpus = { task.attempt <=1 ? '1' : '4' }
    
    input:
		tuple val(chrom), val(filename), val(experiment_name), file('contactCounts.gz'),file('biasFile.gz'), file('mappabilityFile.gz') from contacts_bias_ch.join(map_fithic_ch, by:[0,1,2])

	output:
		tuple val(chrom), val(filename), val(experiment_name), file("*significances.txt.gz") into merge_ch
		file("${filename}*")
	
	conda '/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/conda_envs/fithic'
	
	script:
	"""		
	fithic -i contactCounts.gz -f mappabilityFile.gz -t biasFile.gz	-r ${params.resolution} -o . -l $filename -U -1 -tL 0.5 -tU 6 -v
	"""
}

process mergeContacts {
    tag "${filename}_merge"
    publishDir "${params.outDir}/merged", mode: "copy", overwrite: true

    input:
		tuple val(chrom), val(filename), val(experiment_name), file('significances.txt.gz') from merge_ch
	output:
		tuple val("${experiment_name}_allchr_${params.resolution}_merged.gz"), file("${filename}_merged.gz") into merge_out_ch
	
	conda '/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/conda_envs/fithic'
	
	shell:
	'''
	zcat significances.txt.gz| awk -v q=!{params.fdr} '{if(NR>1 && $7<=q) {print $0}}' | gzip > !{filename}_fithic_subset.gz
	
	#NOTE calls an updated version of CombineNearbyInteraction (CombineNearbyInteraction_PL.py)
	CombineNearbyInteraction_modified.py -i !{filename}_fithic_subset.gz -H 0 -r !{params.resolution} -o !{filename}_merged.gz -k !{chrom}
	
	'''
}

//collect all the chromosomes from each experiment, and merge them
merge_out_ch
  .collectFile(storeDir:"${params.outDir}/merged")
  .set{merge_adjust_ch}
  
process makeWashU {
	tag "${filename}_washu"
    publishDir "${params.outDir}/washu", mode: "copy", overwrite: true

	input:
		file collect_merged from merge_adjust_ch
	output:
		file("${filename}*")
	
	module 'HTSlib/1.11'
	
	shell:
	filename=collect_merged.simpleName
	'''
	#shift positions to bin start, rather than midpoint
	#convert to bed (longrange) format for washu
	#workaround for bug in nextflow where can't edit zipped files, so header line will appear multiple times
	zcat !{collect_merged} | awk -v q=!{params.washUFilter} -v d=!{params.distFilter} '{OFS="\t"} {if($2!="mid1" && $7<q && $4-$2>d) {print $0}}' | awk -v r=!{params.resolution} '{ if ($7==0) {print $1"\t"$14"\t"$15"\t"$3":"$16"-"$17",400"} else {print $1"\t"$14"\t"$15"\t"$3":"$16"-"$17","(-log($7)/log(10))} }' | sort -k1,1 -k2,2n > !{filename}_washu.bed
	bgzip !{filename}_washu.bed
	tabix -f -p bed !{filename}_washu.bed.gz
	'''
}

