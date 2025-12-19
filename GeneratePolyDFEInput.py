import dadi, os
import numpy as np

### First run Rscript_GenerateGeneBinSFS to generate files with the full SFS: each line is a group
### Then use this script to produce PolyDFE input files per group- Stats should match the names of the files,
### which should all have the same filepath, as specified in synfile and nonsynfile.

# for example:
Stats = ['SIFTpolcounts', 'fpkmpolcounts', 'lengthpolcounts', 'NetworkConnectivitypolcounts']


group_number = '15'

for Stat in Stats:

	synfile = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped'+group_number+'_'+Stat+'_synonymous.csv'
	nonsynfile = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/Grouped'+group_number+'_'+Stat+'_nonsynonymous.csv'


	outfile_dir = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome_Clean/PolyDFE'+group_number+'/'
	outfile_start = synfile.split('Grouped'+group_number+'_')[1].split('_')[0]
	outfile_start = outfile_start.split('polcounts')[0]+group_number+'polcounts'
	print(outfile_start)


	###for input files generated from binning, generate a downsampled fs of sample size ns
	###Return statistics for analysis per group: grouper, syn_fs, Ls, Ln
	def fs_generate(inputfile, ns):
		reslist = []
		for line in inputfile.readlines():	
			Ln = line.split(',')[-1]
			grouper = line.split(',')[0]
			SFS = line.split(',')[1:-1]

			SFS.append(0)
			SFS.insert(0,0)
			SFS = np.array(SFS)
		
			syn_fs = dadi.Spectrum(SFS)
			syn_fs_50 = dadi.Spectrum_mod.Spectrum.project(syn_fs, [ns])
			reslist.append([grouper, syn_fs_50, Ln])
		return(reslist)


	###for input files generated from binning, generate an fs.
	###Return statistics for analysis per group: grouper, syn_fs, Ls, Ln
	def fs_generate_full(inputfile):
		reslist = []
		for line in inputfile.readlines():	
			Ln = line.split(',')[-1]
			grouper = line.split(',')[0]
			SFS = line.split(',')[1:-1]
			ns = len(SFS)		
			SFS.append(0)
			SFS.insert(0,0)
			SFS = np.array(SFS)
		
			syn_fs = dadi.Spectrum(SFS)
			reslist.append([grouper, syn_fs, Ln])
		return(reslist)



	with open(synfile, 'r') as synfile, open(nonsynfile, 'r') as nonsynfile:
		reslist_syn = fs_generate_full(synfile)
		reslist_nonsyn = fs_generate_full(nonsynfile)
		reslist_total = zip(reslist_syn, reslist_nonsyn)
		count = 1
		for x in reslist_total:
			SFS_len = str(len(x[0][1])-1)
	# 		with open('/Users/jennyjames/Desktop/ArabidopsisCrossGenome_Clean/SIFTGroup-polcounts'+str(count)+'_'+str(x[0][0])+'_PolyDFE.input', 'a') as resfile:
			with open(outfile_dir+outfile_start+str(count)+'_'+str(x[0][0])+'_PolyDFE.input', 'a') as resfile:
				resfile.write('1 1 '+SFS_len+'\n')
				resfile.write('\t'.join([str(y) for y in x[0][1] if y != 'masked']) + '\t' + str(x[0][-1]))
				resfile.write('\t'.join([str(y) for y in x[1][1] if y != 'masked']) + '\t' + str(x[1][-1]))
				count = count+1
















