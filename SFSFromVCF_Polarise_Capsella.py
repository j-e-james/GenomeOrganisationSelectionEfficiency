import gzip, os, re, numpy, time, random
from collections import defaultdict

### This is a specific script for tackling the Capsella data used in this analysis- ancestral state information
### was inferred from the C. orientalis samples also in the file.

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


def testprint(string):
	if test ==True:
	   print(string)
	else:
		pass

test = True

fold = '4fold'
   
vcf_poly = '/Volumes/MyPassport/Capsella_AlexM/JGIfilesCaprub1_0/Phytozome/PhytozomeV10/Crubella/annotation/corientalis_grandiflora.SNPS.exonic.PacIDs.'+fold+'.vcf.gz'
res_filename = vcf_poly[:-6]+'SFS.csv'

### we can get our monomorphic 0 and 4 fold sites, per pacid protein, from this file (to an approximation):
PacRes = '/Volumes/MyPassport/Capsella_AlexM/JGIfilesCaprub1_0/Phytozome/PhytozomeV10/Crubella/annotation/Crubella_183_v1.0_degeneracy_annotation_PacIDs.txt'

Total_fold_counts = defaultdict(list)
with open(PacRes, 'r') as MonoFile:
	for line in MonoFile:
		Total_fold_counts[line.split('\t')[0].split('pacid=')[1]].append(line.split('\t')[-1][:-1])
print('Loaded annotated file')
		
Genes_dict = defaultdict(list)

### first 16 samples in the VCF are orientalis, the remainder are grandiflora
### we will use this outgroup to assign ancestral states, and just assume parsimony.
orientalis_samples = 16
grandiflora_samples = 50

### Subsample to 50, polydfe can't take that many sampled individuals, and it will match A thaliana...
### we will randomly sample down to this size, to account for missingness due to lack of sequencing information 
subsample_number = 50
subsample = True

with gzip.open(vcf_poly, 'rt') as vcffile:
	for line in vcffile:
		if '#' not in line:
# 			print(line)
# 			if len(line.split('\t')) != 75:
# 				print(line)
			
			#pacid
			pacid= line.split('pacid=')[1].split('\t')[0]
			print(pacid)
			
			### orientalis
			orientalis = ' '.join(line.split('\t')[8:25])
			### we don't care whether data is phased or not
			orientalis = re.findall(' 1[/|\|]0| 0[/|\|]1| 1[/|\|]1| 1[/|\|]0| 0[/|\|]0', orientalis)
			print(orientalis)
			allele_count_orientalis = ''.join(orientalis).count('1')	
			print(allele_count_orientalis)
			Anc = 'not known'
			if allele_count_orientalis == 16*2:
				Anc = 'Alt'
			elif allele_count_orientalis == 0:
				Anc = 'Ref'
				
			if Anc != 'not known':	
				
				###grandiflora
				grandiflora = ' '.join(line.split('\t')[24:])
				### we don't care whether data is phased or not
				grandiflora = re.findall(' 1[/|\|]0| 0[/|\|]1| 1[/|\|]1| 1[/|\|]0| 0[/|\|]0', grandiflora)
	# 			print(grandiflora)
				allele_count_grandiflora = ''.join(grandiflora).count('1')	
				if Anc == 'Alt':
					print('reverse count')
					allele_count_grandiflora = grandiflora_samples*2 - allele_count_grandiflora				
				print(allele_count_grandiflora)	

				if allele_count_grandiflora != 0 and allele_count_grandiflora != grandiflora_samples*2:
					if subsample = True:
				
						"""sample down to subsample_number, and account for missingness. Any site not sequenced in at least subsample_number individuals is removed.
						   It is possible for this subsampling to reduce the number of polymorphic sites observed to 0, so they will not contribute to the SFS."""
	
						allele_count_grandiflora = numpy.random.hypergeometric(allele_count_grandiflora, grandiflora_samples*2-allele_count_grandiflora, subsample_number)
						if allele_count_grandiflora != 0: 
						
							Genes_dict[pacid].append(allele_count_grandiflora)
							
					else:									
						Genes_dict[pacid].append(allele_count_grandiflora)
					
					print(allele_count_grandiflora)	

										
			else:
				print('line skipped')
				
				
for k, v in Genes_dict.items():
	print(k)
	L = Total_fold_counts[k].count(fold)
	print(v)
	SFS = SFS_convert(v, subsample_number)
	print(k+','+','.join(str(x) for x in SFS)+','+str(L))

	with open(res_filename, 'a') as resfile:
		resfile.write(k+','+','.join(str(x) for x in SFS)+','+str(L)+'\n')



