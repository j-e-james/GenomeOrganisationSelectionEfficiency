import gzip, os, re, numpy, time, random
from collections import defaultdict


def testprint(string):
	if test ==True:
	   print(string)
	else:
		pass
   
vcf_poly = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome/1001genomes_snp-short-indel_with_tair10_only_ACGTN.0fold.Poly.Bi.alignmentpolarized.vcf.gz'
vcf_mono = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome/1001genomes_snp-short-indel_with_tair10_only_ACGTN.0fold.Mono.Biallelic.vcf.gz'
genefile = '/Users/jennyjames/Dropbox/AssignAncestral/TAIR10_GFF3_genes.gff.txt.gz'

"""This test value is given to VCFByGene- it throttles the number of lines read in from the VCF
   when splitting based on gene"""
test = False

### will we take a haploid sample
haploid = True

### we will randomly sample down to this size, to account for missingness due to lack of sequencing information 
subsample_number = 50

if haploid:
	res_filename = vcf_poly[:-6]+str(subsample_number)+'.SFS.haploid.csv'
else:
	res_filename = vcf_poly[:-6]+str(subsample_number)+'.SFS.csv'


def SFS_convert(pop_mutations, n):
	"""converts mutation output into SFS"""
	SFS = {} # full SFS with empty sites
	for x in range(1, n):
		SFS[x] = []
	for mut in pop_mutations:
		SFS[int(mut)].append(mut)
	SFS_counts = []
	for k, v in SFS.items():
		SFS_counts.append(len(v))
	return SFS_counts


def GeneSFSProcess(VCFbyGene, gene):
	group_dict = defaultdict(list)
	"""processes a list of VCF lines, that are all the entries for a single gene"""
	for line in VCFbyGene:
		Ref, Alt = line.split('\t')[3], line.split('\t')[4]
		try:
			AncState = line.split('AA=')[1].split(';')[0]
			
		except:
			print('this is our monomorphic file, gene length = '+str(len(VCFbyGene)))
			group_dict[gene].append(len(VCFbyGene))
			"""And we don't need to continue any line by line analyses'"""
			break
			
		testprint(Ref+' '+Alt+' '+AncState)

		### we don't care whether data is phased or not
		polymorphic = re.findall('1[/|\|]0|0[/|\|]1|1[/|\|]1', line)
		not_polymorphic = re.findall('0[/|\|]0', line)
	
		if haploid:
			"""" choose one chromosome randomly per individual. This is relevant here due to very high levels of inbreeding, all are homozygous, i.e. allele_count % 2 == 0 """
		
			### we don't care whether data is phased or not		
			poly_haploid_sample = [re.split('/|\|', x)[random.getrandbits(1)] for x in polymorphic]

			allele_count = ''.join(poly_haploid_sample).count('1')
			polymorphic = allele_count
			not_polymorphic = len(not_polymorphic)

		else:
			### don't sample per individual
			allele_count = ''.join(polymorphic).count('1')	
			polymorphic = allele_count
			not_polymorphic = len(not_polymorphic)*2				

		"""sample down to subsample_number, and account for missingness. Any site not sequenced in at least subsample_number individuals is removed.
		   It is possible for this subsampling to reduce the number of polymorphic sites observed to 0, so they will not contribute to the SFS."""
		if not_polymorphic+polymorphic > subsample_number:
			x = numpy.random.hypergeometric(polymorphic, not_polymorphic, subsample_number)
			polymorphic, not_polymorphic = x, subsample_number-x

			###if resampling hasn't reduced the number of polymorphic sites to 0
			if polymorphic != 0 and not_polymorphic != 0:
				if AncState == Alt:
					derived_state_counts, anc_state_counts = not_polymorphic, polymorphic
				elif AncState == Ref:
					derived_state_counts, anc_state_counts = polymorphic, not_polymorphic
				else:
					print("Unidentified Anc state, (Ref, Alt, AA): "+Ref, Alt, AncState)
					break
				testprint((anc_state_counts, derived_state_counts))
				group_dict[gene].append(derived_state_counts)
			else:
				testprint("Resampling has gotten rid of the site")
				testprint((polymorphic, not_polymorphic))
		else:
			testprint("Sampled in less than half of the individuals")
			testprint(not_polymorphic+polymorphic)
			
	return(group_dict)



def AnalyseVCFByGene(vcf_file, test):

	Genes_dict = defaultdict(list)

	"""Sort VCF by gene"""
	Gene_list = []
	VCF_lines = []
	
	count = 0
	with gzip.open(vcf_file, 'rt') as vcffile:
		for line in vcffile:
			if '#' not in line:
				count = count+1
				Chr, Pos = line.split('\t')[0], int(line.split('\t')[1])
				
				for key in Gene_dict[Chr]:
					if key[0] < Pos < key[1]:
						Gene = Gene_dict[Chr][key]
# 						print(Gene)

						if len(Gene_list) == 0:
							Gene_list.append(Gene)
							VCF_lines.append(line)
							
						elif Gene not in Gene_list:
						
							"""Here we process the previous gene"""							
							print(Gene_list[-1])
							per_gene_dict = GeneSFSProcess(VCF_lines, Gene_list[-1])
							Genes_dict.update(per_gene_dict)

							"""Then we continue, adding our new gene to the list"""
							Gene_list.append(Gene)
							VCF_lines = []
							VCF_lines.append(line)
							
						else:
							VCF_lines.append(line)		

			if test == True:
				if count > 500:
					break
					
		"""Here we process the final gene"""		
		print(Gene_list[-1])
		per_gene_dict = GeneSFSProcess(VCF_lines, Gene_list[-1])
		Genes_dict.update(per_gene_dict)
		
	return(Genes_dict)
 



Gene_dict = defaultdict(dict)
"""GFF3 parser- generate nested dictionary of chromosomes and gene locations"""
with gzip.open(genefile, 'rt') as GFF3_file:
	for line in GFF3_file:
		if line.startswith("#"):
			continue
		elif 'CDS\t' in line:
			try:
				sline = line.split()
				scaf, type, start, end, dir, Parent_info  = sline[0], sline[2], int(sline[3]), int(sline[4]), sline[6], sline[8]
				###Check the use of this: should return the identity of the gene.
				Parent_info = Parent_info.split(";")[0].replace("ID=","").split(",")[0]
				
				Gene_dict[scaf.split('Chr')[1]][start, end] =  (Parent_info.split('Parent=')[1])
				
			except (ValueError, IndexError):
				print("Line conversion failed. Skipping %s.\n" % line)
				continue  



Poly_VCF = AnalyseVCFByGene(vcf_poly, test)
Mono_VCF = AnalyseVCFByGene(vcf_mono, test)

for k, v in Poly_VCF.items():
	print(k)
	print(v)
	SFS = SFS_convert(v, subsample_number)
	testprint(Mono_VCF[k])
	try:
		L = int(Mono_VCF[k][0]) + len(v)
	except:
		"""You won't hit this condition unless testing'"""
		L = 0 +	len(v)
	testprint(L)
	testprint(k+','+','.join(str(x) for x in SFS)+','+str(L))
	with open(res_filename, 'a') as resfile:
		resfile.write(k+','+','.join(str(x) for x in SFS)+','+str(L)+'\n')









	
