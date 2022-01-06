# Index all bam files present in the current directory that start with "filtered_GU" prefix.
ls *.bam | parallel -j 50 samtools index '{}'

# convert all bam to sam in a directory
parallel -j 5 'samtools view -h {} -o {}.sam' ::: *bam

#run mpileup on all bams and save output somewhere else 
parallel -j 5 'samtools mpileup -f /groups/wyattgrp/users/amunzur/mt-ctDNA/references/hg38_mito.fa {} -o /groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/metrics/mpileup/mt_bams_filtered/{}' ::: *.bam
parallel -j 5 'samtools mpileup -f /groups/wyattgrp/users/amunzur/chip_project/references/hg38.fa {} -o /groups/wyattgrp/users/amunzur/chip_project/metrics/mpileup/new_chip_panel_CP/{}' ::: *.bam

# Run the script on all bam files present in the current directories that start with "GU" prefix. 
ls GU*.bam | parallel -j6 -k bash ../process_samples.sh

# List all the bam files in the current directory but ignore them if they start have the string "tmp" in them. 
ls *.bam --ignore="*tmp"

# How to do perl pie on if the string we want to add is a path that contains a slash, this essentially turns bam names into absolute paths.
perl -pi -e 's{filtered_}{/groups/wyattgrp/users/amunzur/chip_project/finland_bams/GU_finland_download/filtered_}g' GU_bam_paths

# move files in the pwd that have a certain prefix to another location 
for file in GU-*; do mv "$file" "/groups/wyattgrp/users/amunzur/chip_project/finland_bams/ctDNA_prognosis_ORIGINAL"; done;

# delete empty lines from a text file with grep 
grep -v '^$' File > newFile

# find the difference in between two files, if neither are sorted: 
comm -2 -3 <(ls | sort) <(awk '{print $2}' ../../data/CHIP_list_of_samples.txt | tail -n +2 | sort)

# start an ssh tunnel to the server to use rmate
ssh -R 52698:localhost:52698 amunzur@lin1

# convert files from DOS-style line endings (CR-LF) to Unix-style line endings
perl -pi -e 's/\r\n/\n/g' input

# remove a string from file names from all files in the current directory, here we remove the string "gz."
for i in *gz.tbi
do
    mv "$i" "`echo $i | sed 's/gz.//'`"
done

# count the number of lines that dont start with #
for file in *vcf
do
	grep "^[^#;]" $file | wc -l 
done

# to copy files over from Finland: 
vpc@tambio.uta.fi:

# how to use pipe to move files: 
ls -thor | grep "Jul 16" | awk '{print $8}' | xargs mv -t ../new_finland_download/
find . -name "*new_finland_download*" -type f | xargs mv -t new_finland_download/

# change extension of many files at once 
for file in *.bam ; do mv $file ${file//bam/mpileup} ; done

# convert VarScan2 output to annovar input 
cat file | awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' | cut -f1,2,3,4,5 | tail -n +2 > output

# collapse strings from a col, useful for grep
cat chip_muts_locations.tsv | cut -f4 | tr '\n' '|'

# after catting many files together remove the repeated header
sed '2,${/^file_name/d;/^---/d}' filename
awk 'FNR==1 { header = $0; print } $0 != header' file > newfile # only keeps the header at the first line of the file

# rename the header of the last col 
for file in *
do
	awk -F, 'NR==1{print "file_name",$0;next}{print FILENAME , $0}' $file > temp
	rm $file
	mv temp $file
done

# perl pie file names to feed into elie's snv calling
ls *.bam | grep WBC | sort | perl -lpe 's/-WBC.*//g'

# make bam list to run in house mut calling pipeline
paste <(ls *.bam | grep WBC | sort | perl -lpe 's/-WBC.*//g') <(ls *.bam | grep -v WBC | sort) <(ls *.bam | grep  WBC | sort) <(find $PWD -name "*.bam" | grep -v WBC | sort) <(find $PWD -name "*.bam" | grep WBC | sort) > tumor_normal_pairs_pipeline.txt

# find the difference between the contents of two files and modify the output with grep
dir1="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures_bowtie/chosen_cfDNA_samples" # chosen 
dir2="/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/figures_bowtie/insert_size_mito_filtered_PNG2" # all 
diff -bur ${dir1} ${dir2} | grep -v "WBC\|Tumor\|CTRL" | cut -d":" -f2 | less -S > samples_removed

# make a bed file by selecting cols using awk and inserting tabs in between
cat edsf_1000012543_hg38_19Jul2021_capture_targets.bed | awk -v OFS='\t' '{print $1,$2,$3}' > capture_targets.bed

# run igv using a batch script
/home/amunzur/IGV_Linux_2.9.4/igv.sh --batch /groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/IGV_batch_scripts/new_chip_panel_ALERT.txt

# run snakemake on the server
snakemake --use-conda --keep-going --cluster-config config/cluster.yaml --cluster 'sbatch -t {cluster.time} --mem={cluster.mem} -c {cluster.cpus} -N {cluster.nodes} -o {cluster.output}' -j 40

# render a Rmd file on terminal
cd /groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/visualization
Rscript -e "rmarkdown::render('interactive_plotting.Rmd')"

# connecting to sftp
sftp gsc-877@sftp.bcgsc.ca