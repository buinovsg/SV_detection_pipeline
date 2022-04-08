#! /usr/bin/env nextflow

params.reference = "/crex/proj/sllstore2017050/nobackup/greta/raw_data/Pvulgaris_442_v2.0.fa"
params.bam_dir = "/crex/proj/sllstore2017050/nobackup/greta/benchmarking/05_pipeline/pipeline_test/bam_dir_test/"
params.chr_list = "/crex/proj/sllstore2017050/nobackup/greta/benchmarking/05_pipeline/pipeline_test/chr_list.txt"
params.vcf_dir = "vcf_output"
params.threads = 8


process runDelly {

	module 'delly'
	module 'bcftools'
        cpus params.threads
	
	publishDir 'params.vcf_dir'
	
	output:
	file "delly.vcf" into delly_vcf_ch
	
	shell:
	'''
	cat !{params.reference} | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t0\t"; } \
	$0 !~ ">" {c+=length($0);} END { print c; }' | tail -n +2 > reference.bed
	grep -wvf !{params.chr_list} reference.bed > exclude_chr_list.bed
	
	export OMP_NUM_THREADS=!{params.threads}
	
	#sample calling
	for bam_file in !{params.bam_dir}*.bam; do
		filename=$(basename $bam_file .bam)
        	for sv_type in {DEL,INS};do
			delly call -t $sv_type -g !{params.reference} -x exclude_chr_list.bed -o "$filename"_"$sv_type".bcf $bam_file
		done
	done
	
	#combined sites bcf
	del_bcf="$(ls *_DEL.bcf)"
	delly merge -t DEL -o sites_DEL.bcf $del_bcf
	ins_bcf="$(ls *_INS.bcf)"
	delly merge -t INS -o sites_INS.bcf $ins_bcf

	#call sv again using sites bcf
	for bam_file in !{params.bam_dir}*.bam; do
        	filename=$(basename $bam_file .bam)
	        for sv_type in {DEL,INS};do
	                delly call -t $sv_type -g !{params.reference} -x exclude_chr_list.bed -v sites_"$sv_type".bcf -o "$filename"_"$sv_type".geno.bcf $bam_file
	        done
	done
	
	#merge sample bcf files
	for sv_type in {DEL,INS};do
	        geno_bcf_files="$(ls *"$sv_type".geno.bcf)"
	        bcftools merge -m id -O b -o merged_"$sv_type".bcf $geno_bcf_files
		bcftools index merged_"$sv_type".bcf
	done
	
	bcftools concat -a merged_*.bcf -o delly.vcf
	'''	
}

process runDysgu {

        cpus params.threads

        publishDir 'params.vcf_dir'

        output:
        file "dysgu.vcf" into dysgu_vcf_ch

        shell:
        '''
        pyenv global 3.8.1

        for bam_file in !{params.bam_dir}*.bam; do
                filename=$(basename $bam_file .bam)
                dysgu run -p!{params.threads} !{params.reference} temp_dir $bam_file > "$filename".vcf
                rm -rf temp_dir
                { grep '#' "$filename".vcf; grep  $'\tPASS\t' "$filename".vcf; } > "$filename"_filter.vcf
        done

        dysgu merge *_filter.vcf > dysgu.vcf
        '''

}

process runManta {

        module 'manta'
        module 'bcftools'
        cpus params.threads

        publishDir 'params.vcf_dir'

        output:
        file "manta.vcf" into manta_vcf_ch

        shell:
        '''
        bams=$(ls "!{params.bam_dir}"*.bam | awk '{print "--bam "$0}' | tr '\n' ' ')
        chrs=$(awk '{print "--region "$0}' !{params.chr_list} | tr '\n' ' ')

        configManta.py --referenceFasta !{params.reference} $bams $chrs
        cd MantaWorkflow
        python runWorkflow.py -m local -j !{params.threads}

        #filter for sites with FILTER set as  PASS
        bcftools view -i "%FILTER='PASS'" results/variants/diploidSV.vcf.gz > ../manta.vcf

        '''
}



delly_vcf_ch
        .concat( dysgu_vcf_ch,manta_vcf_ch )
	.set { vcf_ch }
vcf_ch.view()
