#!/usr/bin/zsh
# apply trimmomatic and salmon quant
# deal with multiple samples at the same time

### fastqc (using R to call fastqc)


### trimmomatic


#<<'COMMENTS'

trim_SE_folder () {
	# a function to trim FASTQ reads using trimmomatic
	# applicable to a folder containing SE FASTQs
	# usage: trim_SE_folder "path_to_a_folder_containing_se_fastq_files"
	
	echo $1

	if [[ -e $1/trimmed/ ]]; then {
		rm -r $1/trimmed/
	} fi

	mkdir $1/trimmed/

	for fastq in `ls $1/|grep ".fastq.gz"`; do {
		fastq_itself=$fastq[0,-10]
		nohup java -jar $SOFT/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 8 \
		-trimlog $1/trimmed/$fastq_itself.log $1/$fastq $1/trimmed/$fastq_itself.trimmed.fastq.gz \
		ILLUMINACLIP:$SOFT/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 \
		SLIDINGWINDOW:4:20 MINLEN:25 \
		&
		if [[ `ps -aux | grep "vrw936" | grep "trimmomatic" | wc -l` > 6 ]]; then {
			wait
		} fi
	} done
	
	wait
	rm $1/trimmed/*.log
}

#COMMENTS




#<<'COMMENTS'

trim_PE_folder () {
	# a function to trim FASTQ reads using trimmomatic
	# applicable to a folder containing PE FASTQs
	# usage: trim_PE_folder "path_to_a_folder_containing_pe_fastq_files"
	
	echo $1

	if [[ -e $1/trimmed/ ]]; then {
		rm -r $1/trimmed/
	} fi

	mkdir $1/trimmed/

	fastq_itself="init"
	for fastq in `ls $1/|grep ".fastq.gz"`; do {
		if [[ "$fastq_itself" != "$fastq[0,-12]" ]]; then {
			fastq_itself=$fastq[0,-12]
			nohup java -jar $SOFT/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 \
			-trimlog $1/trimmed/$fastq_itself.log $1/${fastq_itself}_1.fastq.gz $1/${fastq_itself}_2.fastq.gz $1/trimmed/${fastq_itself}_1.paired.fastq.gz $1/trimmed/${fastq_itself}_1.unpaired.fastq.gz $1/trimmed/${fastq_itself}_2.paired.fastq.gz $1/trimmed/${fastq_itself}_2.unpaired.fastq.gz \
			ILLUMINACLIP:$SOFT/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 \
			SLIDINGWINDOW:4:20 MINLEN:25 \
			&
			if [[ `ps -aux | grep "vrw936" | grep "trimmomatic" | wc -l` > 6 ]]; then {
				wait
			} fi
		} fi
	} done

	wait
	rm $1/trimmed/*.unpaired.fastq.gz
	rm $1/trimmed/*.log
}

#COMMENTS


# GSE125583
# trim_SE_folder "$DATA/ToolTesting/trimmomatic_test/GSE125583";
trim_SE_folder $DATA/scratch/rna_seq_for_mrna/GSE125583;

# GSE125050
# trim_SE_folder "$DATA/ToolTesting/trimmomatic_test/GSE125050";
trim_SE_folder $DATA/scratch/rna_seq_for_mrna/GSE125050;

# GSE110731
# trim_SE_folder "$DATA/ToolTesting/trimmomatic_test/GSE110731";
trim_SE_folder $DATA/scratch/rna_seq_for_mrna/GSE110731;

# GSE104704
# trim_SE_folder "$DATA/ToolTesting/trimmomatic_test/GSE104704";
trim_SE_folder $DATA/scratch/rna_seq_for_mrna/GSE104704;

# GSE95587
# trim_PE_folder "$DATA/ToolTesting/trimmomatic_test/GSE95587";
trim_PE_folder $DATA/scratch/rna_seq_for_mrna/GSE95587;

# GSE53697
# trim_PE_folder "$DATA/ToolTesting/trimmomatic_test/GSE53697";
trim_PE_folder $DATA/scratch/rna_seq_for_mrna/GSE53697;


### salmon quant



#<<'COMMENTS'

salmon_SE_folder () {
	# a function to quantify the expression level of transcripts using salmon quant
	# applicable to a folder containing SE FASTQs
	# usage: salmon_SE_folder "path_to_a_folder_containing_se_fastq_files"
	
	echo $1

	if [[ -e $1/salmon/ ]]; then {
		rm -r $1/salmon/
	} fi

	mkdir $1/salmon/

	for fastq in `ls $1/|grep ".fastq.gz"`; do {
		fastq_itself=$fastq[0,-18]
		nohup salmon quant -i /binf-isilon/alab/students/vrw936/scratch/reference_ome_and_anno \
		-l A -r $1/$fastq \
		--validateMappings --gcBias \
       		-o $1/salmon/$fastq_itself \
		&
		if [[ `ps -aux | grep "vrw936" | grep "salmon quant" | wc -l` > 2 ]]; then {
			wait
		} fi
	} done

	wait  # make sure when the nohup shell_script.sh & is done, all salmon quants are actually done
}

#COMMENTS



#<<'COMMENTS'

salmon_PE_folder () {
	# a function to quantify the expression level of transcripts using salmon quant
	# applicable to a folder containing PE FASTQs
	# usage: salmon_PE_folder "path_to_a_folder_containing_pe_fastq_files"
	
	echo $1

	if [[ -e $1/salmon/ ]]; then {
		rm -r $1/salmon/
	} fi

	mkdir $1/salmon/
	fastq_itself="init"

	for fastq in `ls $1/|grep ".fastq.gz"`; do {
		if [[ "$fastq_itself" != "$fastq[0,-19]" ]]; then {
			fastq_itself=$fastq[0,-19]
			nohup salmon quant -i /binf-isilon/alab/students/vrw936/scratch/reference_ome_and_anno \
			-l A -1 $1/${fastq_itself}_1.paired.fastq.gz \
			-2 $1/${fastq_itself}_2.paired.fastq.gz \
			--validateMappings --gcBias \
			-o $1/salmon/$fastq_itself \
			&
			if [[ `ps -aux | grep "vrw936" | grep "salmon quant" | wc -l` > 2 ]]; then {
				wait
			} fi
		} fi
	} done

	wait
}

#COMMENTS

# GSE125583
# salmon_SE_folder $DATA/ToolTesting/trimmomatic_test/GSE125583/trimmed
salmon_SE_folder $DATA/scratch/rna_seq_for_mrna/GSE125583/trimmed;

# GSE125050
# salmon_SE_folder $DATA/ToolTesting/trimmomatic_test/GSE125050/trimmed;
salmon_SE_folder $DATA/scratch/rna_seq_for_mrna/GSE125050/trimmed;

# GSE110731
# salmon_SE_folder $DATA/ToolTesting/trimmomatic_test/GSE110731/trimmed;
salmon_SE_folder $DATA/scratch/rna_seq_for_mrna/GSE110731/trimmed;

# GSE104704
# salmon_SE_folder $DATA/ToolTesting/trimmomatic_test/GSE104704/trimmed;
salmon_SE_folder $DATA/scratch/rna_seq_for_mrna/GSE104704/trimmed;

# GSE95587
# salmon_PE_folder $DATA/ToolTesting/trimmomatic_test/GSE95587/trimmed
salmon_PE_folder $DATA/scratch/rna_seq_for_mrna/GSE95587/trimmed;

# GSE53697
# salmon_PE_folder $DATA/ToolTesting/trimmomatic_test/GSE53697/trimmed;
salmon_PE_folder $DATA/scratch/rna_seq_for_mrna/GSE53697/trimmed;


################################ AMP-AD datasets ##################################


##### rename files

# ROSMAP
ls $DATA/scratch/rna_seq_for_mrna/ROSMAP/fastq/*.fastq.gz > $DATA/scratch/rna_seq_for_mrna/ROSMAP/fastq/former_file_names.txt

for fastq in `ls $DATA/scratch/rna_seq_for_mrna/ROSMAP/fastq |grep ".r[1|2].fastq.gz"`; do { 
	mv $DATA/scratch/rna_seq_for_mrna/ROSMAP/fastq/$fastq $DATA/scratch/rna_seq_for_mrna/ROSMAP/fastq/${fastq[0,-13]}_${fastq[-10]}.fastq.gz 
} done

# MSBB
#ls $DATA/scratch/rna_seq_for_mrna/MSBB/fastq/*.fastq.gz > $DATA/scratch/rna_seq_for_mrna/MSBB/fastq/former_file_names.txt

for fastq in `ls $DATA/scratch/rna_seq_for_mrna/MSBB/fastq |grep ".accepted_hits.sort.coord.combined.fastq.gz"`; do { 
	mv $DATA/scratch/rna_seq_for_mrna/MSBB/fastq/$fastq $DATA/scratch/rna_seq_for_mrna/MSBB/fastq/${fastq[0,-44]}.fastq.gz 
} done

# MayoRNAseq

# CBE
#ls $DATA/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_CBE_sample_FASTQs/*.fastq.gz > $DATA/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_CBE_sample_FASTQs/former_file_names.txt

for fastq in `ls $DATA/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_CBE_sample_FASTQs |grep ".r[1|2].fastq.gz"`; do {
	mv $DATA/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_CBE_sample_FASTQs/$fastq $DATA/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_CBE_sample_FASTQs/${fastq[0,-35]}_${fastq[-10]}.fastq.gz
} done


# TCX

#ls $DATA/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_TCX_sample_FASTQs/*.fastq.gz > $DATA/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_TCX_sample_FASTQs/former_file_names.txt

for fastq in `ls $DATA/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_TCX_sample_FASTQs |grep ".r[1|2].fastq.gz"`; do {
	mv $DATA/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_TCX_sample_FASTQs/$fastq $DATA/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_TCX_sample_FASTQs/${fastq[0,-13]}_${fastq[-10]}.fastq.gz
} done





##### trim, shuffle and quant

shuffle_PE(){
	# a function to randomise to order of reads in the given PE fastq files
	# usage example: shuffle_PE $DATA/scratch/rna_seq_for_mrna/ROSMAP/fastq/trimmed

	echo $1

	if [[ -e $1/shuffled/ ]]; then {
		rm -r $1/shuffled/
	} fi

	mkdir $1/shuffled/
	fastq_itself="init"

	for fastq in `ls $1/|grep ".fastq.gz"`; do {
		if [[ "$fastq_itself" != "$fastq[0,-19]" ]]; then {
			fastq_itself=$fastq[0,-19]
			nohup shuffle.sh in1=$1/${fastq_itself}_1.paired.fastq.gz in2=$1/${fastq_itself}_2.paired.fastq.gz out1=$1/shuffled/${fastq_itself}_1.paired.fastq.gz out2=$1/shuffled/${fastq_itself}_2.paired.fastq.gz &
			if [[ `ps -aux | grep "vrw936" | grep "shuffle.sh" | wc -l` > 1 ]]; then {
				wait
			} fi
		} fi
	} done

	wait
}



shuffle_SE(){
	# a function to randomise to order of reads in the given PE fastq files
	# usage example: shuffle_PE $DATA/scratch/rna_seq_for_mrna/ROSMAP/fastq/trimmed
	echo $1

	if [[ -e $1/shuffled/ ]]; then {
		rm -r $1/shuffled/
	} fi

	mkdir $1/shuffled/

	for fastq in `ls $1/|grep ".fastq.gz"`; do {
		nohup shuffle.sh in=$1/$fastq out=$1/shuffled/$fastq &
		if [[ `ps -aux | grep "vrw936" | grep "shuffle.sh" | wc -l` > 1 ]]; then {
			wait
		} fi
	} done

	wait
}


# ROSMAP
trim_PE_folder $DATA/scratch/rna_seq_for_mrna/ROSMAP/fastq;
shuffle_PE $DATA/scratch/rna_seq_for_mrna/ROSMAP/fastq/trimmed;
salmon_PE_folder $DATA/scratch/rna_seq_for_mrna/ROSMAP/fastq/trimmed/shuffled;

# MSBB
trim_SE_folder $DATA/scratch/rna_seq_for_mrna/MSBB/fastq;
shuffle_SE $DATA/scratch/rna_seq_for_mrna/MSBB/fastq/trimmed;
salmon_SE_folder $DATA/scratch/rna_seq_for_mrna/MSBB/fastq/trimmed/shuffled;


# MayoRNAseq
# CBE
trim_PE_folder $DATA/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_CBE_sample_FASTQs;
shuffle_PE $DATA/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_CBE_sample_FASTQs/trimmed;
salmon_PE_folder $DATA/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_CBE_sample_FASTQs/trimmed/shuffled;

# TCX
trim_PE_folder $DATA/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_TCX_sample_FASTQs;
shuffle_PE $DATA/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_TCX_sample_FASTQs/trimmed;
salmon_PE_folder $DATA/scratch/rna_seq_for_mrna/MayoRNAseq/fastq/Mayo_TCX_sample_FASTQs/trimmed/shuffled;

