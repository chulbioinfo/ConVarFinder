import os, glob, math, copy


def Parse_TAAS_file(file_name):

	gene_taas_dic = {}

	fin = open(file_name,'r')
	fin.readline()
	for line in fin:
		gene,n_pos,a_pos = line.strip().split('\t')
		gene_taas_dic.setdefault(gene,[])
		gene_taas_dic[gene].append([int(n_pos),int(a_pos)])
	fin.close()

	return gene_taas_dic
	


def Parse_BLAST_file(file_name):

	gene_align_dic = {}

	fin = open(file_name,'r')
	for line in fin:

		qid = line.strip().split('\t')[0]
		tid = line.strip().split('\t')[1]
		qst = int(line.strip().split('\t')[7])
		qend = int(line.strip().split('\t')[8])
		tst = int(line.strip().split('\t')[9])
		tend = int(line.strip().split('\t')[10])
		
		gene_align_dic.setdefault(qid,[])
		gene_align_dic[qid].append([qst,qend,tid,tst,tend])

	fin.close()

	for qid in gene_align_dic.keys():
		gene_align_dic[qid] = sorted(gene_align_dic[qid], key=lambda x: (x[0],x[1]))

	return gene_align_dic



def Parse_MSA(dir_name, species_name):

	sequence_dic = {}

	fasta_list = glob.glob(dir_name+'/*.fas')
	for f in fasta_list:

		gene_name = f.split('/')[-1].split('.')[0]
		flag = 0
		fin = open(f,'r')
		for line in fin:

			if '>'+species_name in line:
				flag = 1
				sequence_dic.setdefault(gene_name,'')
			elif '>' in line[0]:
				flag = 0

			elif flag == 1:
				sequence_dic[gene_name]+= line.strip()
		fin.close()

	return sequence_dic 



def Convert_taas_to_gapless_pos(taas_dic,seq_dic):

	New_taas_dic = {}

	for gene in sorted(taas_dic.keys()):

		if gene in seq_dic.keys():
			New_taas_dic.setdefault(gene,[])

			for taas in taas_dic[gene]:
				
				if not seq_dic[gene][(taas[1]-1)*3:taas[1]*3] == '---':

					gapless_n_pos = len(seq_dic[gene][:(taas[1]-1)*3+1].replace('-',''))
					gapless_a_pos = int(math.ceil(float(len(seq_dic[gene][:(taas[1]-1)*3+1].replace('-','')))/3.0))

					New_taas_dic[gene].append([taas[0],taas[1],gapless_n_pos,gapless_a_pos])

	return New_taas_dic



def Convert_taas_to_genome_pos(gapless_taas_dic,blast_f):

	New_taas_dic = {}

	blast_dic = Parse_BLAST_file(blast_f)
	temp_dic = copy.deepcopy(gapless_taas_dic)

	for gene in sorted(temp_dic.keys()):
		
		New_taas_dic.setdefault(gene,[])

		if gene in blast_dic.keys():

			for blast in blast_dic[gene]:

				for taas in temp_dic[gene]:

					gap_npos,gap_apos,npos,apos = taas
					qst,qend,chrom,gst,gend = blast
					new_pos = 'na'

					if qst <= npos <= qend:

						if gend - gst >= 0:
							strand = '+'
							new_pos = gst+(npos-qst)
						elif gend - gst < 0:
							strand = '-'
							new_pos = gst-(npos-qst)
						New_taas_dic[gene].append([gap_npos,gap_apos,npos,apos,chrom,new_pos,strand])

	
	cnt = 0
	new_total_cnt = 0
	total_cnt = 0
	for gene in sorted(gapless_taas_dic.keys()):
		New_taas_string = []
		for taas in New_taas_dic[gene]:
			if not taas[:-3] in New_taas_string:
				New_taas_string.append(taas[:-3])
				new_total_cnt += 1
		for taas in gapless_taas_dic[gene]:
			total_cnt += 1
			if taas not in New_taas_string:
				cnt += 1

	return New_taas_dic



def Compare_with_vcf(taas_dic,vcf_f,OutFile):

	fout = open(OutFile,'w')
	fout.write('Chromosome\tChromosomal_Position\tSNP\tReference_Allele\t'+ \
			   'Alternative_Allele\tMSA_Position(Nuc)\tMSA_Position(AA)\tVCF_Annotation\n')

	taas_pos_dic = {}
	taas_ann_dic = {}
	for gene in sorted(taas_dic.keys()):
		for taas in taas_dic[gene]:
			gap_npos,gap_apos,npos,apos,chrom,new_pos,strand = taas
			taas_pos_dic.setdefault(chrom,[])
			taas_pos_dic[chrom].append(new_pos)
			taas_ann_dic[(chrom, new_pos)] = [gene,gap_npos,gap_apos,npos,apos,chrom,new_pos,strand]
	for chrom in taas_pos_dic.keys():
		taas_pos_dic[chrom] = sorted(taas_pos_dic[chrom])

	fin = open(vcf_f,'r')
	for line in fin:
		if not '#' in line[0]:
			chrom = line.split('\t')[0]
			if chrom in taas_pos_dic.keys():
				pos = int(line.split('\t')[1])
				ref, alt = line.split('\t')[3:5]
				snp_id = line.split('\t')[2]

				if ';ANN=' in line:
					Anns = line.strip().split('\t')[7].split(';ANN=')[1].split(';')[0].split(',')
				else:
					Anns = []
				VcfGenes = ''
				for Ann in Anns:

					VcfGene = Ann.split('|')[3]
					VarType = Ann.split('|')[1]
					VcfGenes += VcfGene+'_'+VarType+';'

				if pos in taas_pos_dic[chrom]:
					gene,gap_npos,gap_apos,npos,apos,chrom,new_pos,strand = taas_ann_dic[(chrom, pos)]

					if 'missense_variant' in VcfGenes:

						fout.write(chrom+'\t'+str(pos)+'\t'+snp_id+'\t'+ref+'\t'+alt+'\t'+ \
							str(npos)+'\t'+str(apos)+'\t'+VcfGenes+'\n')

	fin.close()
	fout.close()



if __name__ == "__main__":

	## Input Files ##

	TaasFile = 'TAAS_avians.txt'
	Chicken_VcfFile = '../eff_Gallus_gallus.vcf' # Needs to be downloaded from "https://drive.google.com/file/d/1jZMFxDjGmRdN5lRWKwVVLDcGrPgg_GPy/view?usp=sharing"
	ZebraFinch_VcfFile = '../eff_Taeniopygia_guttata.vcf' # Needs to be downloaded from "https://drive.google.com/file/d/1jZMFxDjGmRdN5lRWKwVVLDcGrPgg_GPy/view?usp=sharing"

	Chicken_BlastOut = 'GALGA_from_Avian_MSA_Sequence.blast_out'
	ZebraFinch_BlastOut = 'TAEGU_from_Avian_MSA_Sequence.blast_out'

	#################


	## Parse TAAS ##

	Avian_taas_dic = Parse_TAAS_file(TaasFile)


	## Parse Multiple sequence alignment (MSA) File ##

	G_A_seq_dic = Parse_MSA('Avian_MSA_Sequence','GALGA')
	T_A_seq_dic = Parse_MSA('Avian_MSA_Sequence','TAEGU')

	gapless_G_A_taas_dic = Convert_taas_to_gapless_pos(Avian_taas_dic,G_A_seq_dic)
	gapless_T_A_taas_dic = Convert_taas_to_gapless_pos(Avian_taas_dic,T_A_seq_dic)


	## Convert TAAS MSA position to Chromosomal position ##

	g_G_A_taas_dic = Convert_taas_to_genome_pos(gapless_G_A_taas_dic,Chicken_BlastOut)
	g_T_A_taas_dic = Convert_taas_to_genome_pos(gapless_T_A_taas_dic,ZebraFinch_BlastOut)


	## Check if TAAS overlaps with Polymorphisms ##
	
	Compare_with_vcf(g_G_A_taas_dic,Chicken_VcfFile,'FixedDifference_Chicken.txt')
	Compare_with_vcf(g_T_A_taas_dic,ZebraFinch_VcfFile,'FixedDifference_ZebraFinch.txt')
