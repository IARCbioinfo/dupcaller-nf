#! /usr/bin/env nextflow

// Copyright (C) 2010 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.mem           = 32
params.cpu           = 16
params.ref           = null
params.region        = "chr{1..22} chrX chrY"
params.known_snp     = null
params.noise_mask    = null
params.tn_file       = null
params.output_folder = "dupcaller_results"

params.help = null

log.info "" 
log.info "--------------------------------------------------------"
log.info "  dupcaller-nf 1.0.0: Dupcaller pipeline for somatic variant calling with Nextflow "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""


if (params.help) {
log.info '-------------------------------------------------------------'
    log.info ' USAGE  '
    log.info '-------------------------------------------------------------'
    log.info ''
    log.info 'DupCallerCall.py -b ${sample}.bam -f reference.fa -o {output_predix} -p {threads} -n {normal.bam} -g germline.vcf.gz -m noise_mask.bed.gz'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --tn_file            FILE                    input tabulation-separated values file with columns sample (sample name),'
    log.info '                                                 tumor (full path to tumor bam), normal (full path to matched normal bam);'
    log.info '                                                 optionally (for --genotype mode), columns preproc (is the bam RNAseq needing'
    log.info '                                                 preprocessing: yes or no) and vcf (full path to vcf file containing alleles to genotype)'
    log.info '                                                 where each line contains the path of two matched BAM files.'    
    log.info '    --ref                FILE (with indexes)     Reference fasta file.'
    log.info ''
    log.info 'Optional arguments:'
    log.info '    --noise_mask         FILE                    Bed file containing intervals to mask.'
    log.info '    --known_snp          FILE                    VCF file with known variants and frequency (e.g., from gnomad).'
    log.info '    --cpu                INTEGER                 Number of cpu used (default: 4).'
    log.info '    --mem                INTEGER                 Java memory passed to mutect in GB (default: 8).'
    log.info '    --output_folder      FOLDER                  Output folder (default: mutect_results).'
    log.info ''
    exit 0
}else{
    /* Software information */
    log.info "mem                    = ${params.mem}"
    log.info "cpu                    = ${params.cpu}"
    log.info "output_folder          = ${params.output_folder}"
    log.info "noise_mask             = ${params.noise_mask}"
    log.info "known_snp              = ${params.known_snp}"
    log.info "tn_file                = ${params.tn_file}"
    log.info "ref                    = ${params.ref}"
}

// reference file and its indexes
ref = tuple file(params.ref), file(params.ref+'.fai'), file(params.ref+'.gzi'), 
    file( params.ref.replace(".fasta",".dict").replace(".fa",".dict") )

known_snp = params.known_snp ? (tuple file(params.known_snp), file(params.known_snp+".tbi")) : (tuple file("NO_SNP"), file("NO_SNP.tbi"))
noise_mask = params.noise_mask ? (file(params.noise_mask)) : (file("NO_BED"))


// Process to extract chromosome names from the reference index file (.fai)
process get_chromosomes {
    input:
    tuple path(fasta_ref), path(fasta_ref_fai), path(fasta_ref_hp_h5), path(fasta_ref_h5), path(fasta_tn_h5)

    output:
    //val(chroms), emit: chrom_list
    //stdout
    path("chrom_list.txt")

    script:
    """
    cut -f1 ${fasta_ref_fai} | sort | uniq | grep -v -P "alt|random|Un|chrEBV|HLA" > chrom_list.txt
    """

    stub:
    """
    cut -f1 ${fasta_ref_fai} | sort | uniq | grep -v -P "alt|random|Un|chrEBV|HLA" > chrom_list.txt
    """
    
}

//dupCallerCall
process dupCallerCall{

    memory params.mem+'GB'
    cpus params.cpu
    maxForks 15

    publishDir "${params.output_folder}/dupcaller/", mode: 'copy' 

    input:
        tuple val(sample), path(bamT), path(baiT), path(bamN), path(baiN)
        tuple path(fasta_ref), path(fasta_ref_fai), path(fasta_ref_hp_h5), path(fasta_ref_h5), path(fasta_tn_h5)
        tuple path(vcf), path(vcf_tbi)
        path(bed)
        each chromosome


    output:
        //tuple val(sample), path("${sample}*snv.vcf"), path("${sample}*indel.vcf"), emit : vcfs
        tuple val(sample), path("${sample}*/${sample}*_snv.vcf"), path("${sample}*/${sample}*_indel.vcf"), emit : vcfs
        path("${sample}*/*"), emit : report

    shell:
        germline = (vcf.baseName=="NO_VCF") ? "": "-g " + vcf.join(" -g ")
        noise_mask = (bed.baseName=="NO_BED") ? "" : "-m " + bed.join(" -m ")
        """
        DupCaller.py call -tt 30 -b $bamT -n $bamN  -f $fasta_ref -o ${sample}_${chromosome} -p $task.cpus -r $chromosome $germline $noise_mask 
        """

    stub:
        """
        touch ${sample}_${chromosome}/${sample}_${chromosome}.snv.vcf
        touch ${sample}_${chromosome}/${sample}_${chromosome}.indel.vcf
        touch ${sample}_${chromosome}/${sample}_${chromosome}.txt
        touch ${sample}_${chromosome}/${sample}_${chromosome}.png
        """

}



process mergeResults{

    memory params.mem+'GB'
    cpus params.cpu

    publishDir "${params.output_folder}/", mode: 'copy', pattern: '{*vcf.gz}' 

    input:
        tuple val(sample), path(snvs), path(indels)

    output:
        path("${sample}_calls.vcf.gz")

    shell:
        """
        # MERGE VCF FILES
        sed '/^#CHROM/q' `ls -1 *.vcf | head -1` > header.txt

        # Determine sort options based on the availability of version-sort in sort command
        sort_ops=\$(sort --help | grep -q 'version-sort' && echo "-k1,1V" || echo "-k1,1d")

        # Concatenate VCF contents and sort
        grep --no-filename -v '^#' *.vcf | LC_ALL=C sort -T \$PWD -t '	' \$sort_ops -k2,2n >> header.txt

        # Rename the merged file
        mv header.txt ${sample}_calls.vcf

        # add GT tag for annovar
        ${projectDir}/bin/fixDupcallerOutput.sh *_calls.vcf
	
        # compress
        bgzip ${sample}_calls.vcf
        """

    stub:
        """
        touch ${sample}_calls.vcf.gz
        """
}

/***************************************************************************************/
/************************  Workflow : main *********************************************/
/***************************************************************************************/


workflow {

    //load input files
    assert(params.tn_file) : "Error: please provide tn_file"
    pairs = Channel.fromPath(params.tn_file).splitCsv(header: true, sep: '\t', strip: true)
        .map{ row -> 
        assert(row.sample) : "Error: sample tag missing check your input_file"
        assert(row.tumor) : "Error: tumor bam is missing check your input_file"
        assert(row.normal) : "Error: normal bam is missing check your input_file"
        tuple(
            row.sample,
            file(row.tumor), file("${row.tumor}.{bai,crai}")[0],
            file(row.normal), file("${row.normal}.{bai,crai}")[0]
        )}


    //required_extensions = ['0123', 'amb', 'ann', 'bwt.8bit.32', 'pac']
    //index_files = required_extensions.collect { ext -> file("${ref_base}.fasta.${ext}")
    //dup_idx_ch = Channel.of( index_files.collect { file(it) })

    ref = tuple( file(params.ref), 
                file("${params.ref}.fai"), 
                    //file("${params.ref}".replaceAll(/\.fa(sta)?$/, '.dict')),
                file("${params.ref}.hp.h5"),
                file("${params.ref}.ref.h5"),
                file("${params.ref}.tn.h5")
    )

    chromosome = get_chromosomes(ref).splitText().collect{ it.trim() }.view()

    ///dupCallerCall(pairs,ref,params.region,known_snp,noise_mask)
    dupCallerCall(pairs, ref, known_snp, noise_mask, chromosome)
    mergeResults( dupCallerCall.out.vcfs.groupTuple(by: 0) )


}
