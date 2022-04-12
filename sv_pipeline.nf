#! /usr/bin/env nextflow

params.reference = "/crex/proj/sllstore2017050/nobackup/greta/raw_data/Pvulgaris_442_v2.0.fa"
params.bam_dir = "/crex/proj/sllstore2017050/nobackup/greta/benchmarking/05_pipeline/pipeline_test/bam_dir_test/"
params.chr_list = "/crex/proj/sllstore2017050/nobackup/greta/benchmarking/05_pipeline/pipeline_test/chr_list.txt"
params.vcf_dir = "$PWD/vcf_output"
params.threads = 8
params.overlap_size = 0.8

process runDelly {

	module 'delly'
	module 'bcftools'
        cpus params.threads
	
	publishDir params.vcf_dir
	
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
        	for sv_type in {DEL,INS,INV};do
			delly call -t $sv_type -g !{params.reference} -x exclude_chr_list.bed -o "$filename"_"$sv_type".bcf $bam_file
		done
	done
	
	#combined sites bcf
	del_bcf="$(ls *_DEL.bcf)"
	delly merge -t DEL -o sites_DEL.bcf $del_bcf
	ins_bcf="$(ls *_INS.bcf)"
	delly merge -t INS -o sites_INS.bcf $ins_bcf
	inv_bcf="$(ls *_INV.bcf)"
        delly merge -t INV -o sites_INV.bcf $inv_bcf

	#call sv again using sites bcf
	for bam_file in !{params.bam_dir}*.bam; do
        	filename=$(basename $bam_file .bam)
	        for sv_type in {DEL,INS,INV};do
	                delly call -t $sv_type -g !{params.reference} -x exclude_chr_list.bed -v sites_"$sv_type".bcf -o "$filename"_"$sv_type".geno.bcf $bam_file
	        done
	done
	
	#merge sample bcf files
	for sv_type in {DEL,INS,INV};do
	        geno_bcf_files="$(ls *"$sv_type".geno.bcf)"
	        bcftools merge -m id -O b -o merged_"$sv_type".bcf $geno_bcf_files
		bcftools index merged_"$sv_type".bcf
	done
	
	bcftools concat -a merged_*.bcf -o delly.vcf
	'''	
}

process runDysgu {
	
	module 'bcftools'
        cpus params.threads

        publishDir params.vcf_dir

        output:
        file "dysgu.vcf" into dysgu_vcf_ch

        shell:
        '''
        pyenv global 3.8.1

        for bam_file in !{params.bam_dir}*.bam; do
                filename=$(basename $bam_file .bam)
                dysgu run -p!{params.threads} !{params.reference} temp_dir $bam_file > "$filename".vcf
                rm -rf temp_dir
                bcftools view -i "%FILTER='PASS'" "$filename".vcf  > "$filename"_filter.vcf
        done

        dysgu merge *_filter.vcf > dysgu.vcf
        '''

}

process runManta {

        module 'manta'
        module 'bcftools'
        cpus params.threads

        publishDir params.vcf_dir

        output:
        file "manta.vcf" into manta_vcf_ch

        shell:
        '''
        bams=$(ls "!{params.bam_dir}"*.bam | awk '{print "--bam "$0}' | tr '\n' ' ')
        chrs=$(awk '{print "--region "$0}' !{params.chr_list} | tr '\n' ' ')

        configManta.py --referenceFasta !{params.reference} $bams $chrs
        cd MantaWorkflow
        python runWorkflow.py -m local -j !{params.threads}

        bcftools view results/variants/diploidSV.vcf.gz > ../manta.vcf

        '''
}

delly_vcf_ch
        .concat( dysgu_vcf_ch,manta_vcf_ch )
	.set { vcf_ch }
//vcf_ch.view()

process overlap_vcf {

        module 'bcftools'
        module 'BEDTools'

        publishDir baseDir

	input:
	file vcf from vcf_ch.collect()
	
        output:
        file "manta_overlap.vcf"

        shell:
        '''
        #deletion & inversion .bed
        for file in !{params.vcf_dir}/*.vcf; do
                
		filename=$(basename $file .vcf)
                
		bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%FILTER\n' $file | \
                awk '($5 == "DEL" && $6 == "PASS") || ($5 == "DEL" && $6 == ".")' > "$filename"_del.bed
		
		bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%FILTER\n' $file | \
                awk '($5 == "INV" && $6 == "PASS") || ($5 == "INV" && $6 == ".")' > "$filename"_inv.bed
        done

        #insertion .bed
        bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%FILTER\t%INFO/INSLEN\n' !{params.vcf_dir}/delly.vcf | awk '($5 == "INS" && $6 == "PASS") || ($5 == "INS" && $6 == ".")' | awk '{print $0, $2+$7}' | awk -v OFS='\t' '{print $1,$2,$8,$4,$5,$6,$7}' > delly_ins.bed
        bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%FILTER\t%INFO/SVLEN\n' !{params.vcf_dir}/dysgu.vcf | awk '($5 == "INS" && $6 == "PASS") || ($5 == "INS" && $6 == ".")' | awk '{print $0, $2+$7}' | awk -v OFS='\t' '{print $1,$2,$8,$4,$5,$6,$7}' > dysgu_ins.bed
        bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%FILTER\t%INFO/SVLEN\n' !{params.vcf_dir}/manta.vcf | awk '($5 == "INS" && $6 == "PASS") || ($5 == "INS" && $6 == ".")' | awk '{print $0, $2+$7}' | awk -v OFS='\t' '{print $1,$2,$8,$4,$5,$6,$7}' | awk '$7 != "."' > manta_ins.bed


        for sv_type in {del,ins,inv}; do
                bedtools multiinter -i dysgu_"$sv_type".bed delly_"$sv_type".bed manta_"$sv_type".bed > common_overlap_"$sv_type".bed
                bedtools intersect -a common_overlap_"$sv_type".bed -b dysgu_"$sv_type".bed delly_"$sv_type".bed manta_"$sv_type".bed -f !{params.overlap_size} -r -wa -wb > common_overlap_with_ID_"$sv_type".bed
                bedtools groupby -g 1-8 -c 9 -o count_distinct -i common_overlap_with_ID_"$sv_type".bed > distinct_overlap_"$sv_type".bed
                awk '$9==3' distinct_overlap_"$sv_type".bed | awk -v OFS='\t' '{print $1,$2,$3}' > for_matching_"$sv_type".txt
                grep -F -w -f for_matching_"$sv_type".txt common_overlap_with_ID_"$sv_type".bed | awk '$9=="3"' | awk '{print $13}' > ID_list_"$sv_type".txt
        done

        cat ID_list_del.txt ID_list_ins.txt ID_list_inv.txt > ID_list.txt
        bcftools view -i'ID=@ID_list.txt' !{params.vcf_dir}/manta.vcf > manta_overlap.vcf

        '''
}

