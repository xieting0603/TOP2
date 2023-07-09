#!/usr/bin/python3
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; 
#
import math
import straw
import numpy as np
import argparse
import collections
import cooler
parser = argparse.ArgumentParser(description='Slide a square along the diagonal of a log2 ratio map')
parser.add_argument('--hic1',  nargs=1,
                    help='Juicebox hic file 1', required=True)
parser.add_argument('--hic2',  nargs=1,
                    help='Juicebox hic file 2', required=True)
parser.add_argument('--binSize', type=int , nargs=1,
                    help='Binsize of HiC in bp', required=True)
parser.add_argument('--windowSize', type=int , nargs=1,
                    help='Window size in bins (i.e. at windowSize of 10 with 25000 bp binSize means a window size of 250000bp)', required=True)
parser.add_argument('--createHiC', default=False, action='store_true',
                    help='Create hic output files')
parser.add_argument('--exp1',  nargs=1,
                    help='Juicebox hic file 1', required=True)
parser.add_argument('--exp2',  nargs=1,
                    help='Juicebox hic file 2', required=True)

args = parser.parse_args()
windowSize = args.windowSize[0]
exp_file1=args.exp1[0]
exp_file2=args.exp2[0]

hicFile1    = args.hic1[0]
hicFile2    = args.hic2[0]
binSize = args.binSize[0]
binSizekb = int(args.binSize[0])//1000
createHiC   = args.createHiC
prefix=hicFile1.split('_')[0]+'_'+hicFile2.split('_')[0]+'_'+str(binSizekb)+'_w'+str(windowSize)

bedgraphFile   = 'DKOWT_'+str(binSizekb)+'_w'+str(windowSize)+"_g2.bedgraph111"
exp1={}
exp2={}
c=cooler.Cooler(hicFile1+'::/resolutions/'+str(binSize))

chroms=c.chromnames
def obs_exp_matrix(cool_file_mat,exp_file,chrr):
    exp={}
    for i in open(exp_file):
        t=i.strip().split()
        if not t[0]=='chrom' and not t[-1]==np.nan:
            exp[(t[0],int(t[1]))]=t[-1]
    m_exp=np.zeros_like(cool_file_mat)
    for i in range(0,m_exp.shape[0]):
        for j in range(0,m_exp.shape[0]):
            m_exp[i,j]=exp[(chrr,abs(i-j))]
    m_nor=cool_file_mat/m_exp
    m_nor[m_exp == 0 ]==0
    m_nor[np.isnan(cool_file_mat)]==np.nan
    m_nor[np.isnan(m_exp)]==np.nan

    return m_nor

# Sliding Window
with open(bedgraphFile, "w") as f:
	for chrom in chroms:

		c1 = cooler.Cooler(hicFile1+'::/resolutions/'+str(binSize))
		c2 = cooler.Cooler(hicFile2+'::/resolutions/'+str(binSize))
		m_raw1=c1.matrix(balance=True).fetch(chrom)
		m_raw2=c2.matrix(balance=True).fetch(chrom)
		m1=obs_exp_matrix(m_raw1,exp_file1,chrom)
		m2=obs_exp_matrix(m_raw2,exp_file2,chrom)
		assert m1.shape[0] == m1.shape[1]
		assert m1.shape[0] == m2.shape[0]
		assert m1.shape[1] == m2.shape[1]
		#selected overlap non-nan region
		

		m_log2 = np.log2(m1 / m2)
		# Set everywhere 0 where no value is available in either map
		m_log2[m1 == 0 ] = 0
		m_log2[m2 == 0 ] = 0
		m_log2[np.isnan(m1)]==np.nan
		m_log2[np.isnan(m2)]==np.nan
		print("Slide square on chrom %s" % chrom)
		# Sliding window
		for i in range(windowSize, (m_log2.shape[0]-windowSize) ):
			# rows
			s1 = i-windowSize
			e1 = i
			
			# columns
#			s2 = (i+1)
#			e2 = (i+windowSize+1)
			s2 = i-windowSize
			e2 = i

			assert 0 <= s1 and s1 < m_log2.shape[0]
			assert 0 <= e1 and e1 <= m_log2.shape[0]
			assert 0 <= s2 and s2 < m_log2.shape[1]
			assert 0 <= e2 and e2 <= m_log2.shape[1]

			sub = m_log2[s1:e1, s2:e2]
#			sub=sub[sub!=0]
			value = np.nanmean(sub)
			shape=sub.shape[0]
#			print('nan='+str(len(sub[np.isnan(sub)])))
#			print('zero='+str(np.nansum(np.where(sub,0,1))))
			start = (i-windowSize//2)*binSize
			end = start + binSize
			if not value==value or len(sub[np.isnan(sub)])>0.9*shape*shape:
#			if not value==value :
				value=0
			f.write("%s\t%d\t%d\t%s\n" % (chrom, start, end, value))
		if createHiC == True:
			# Write sparse triplet
			with open((prefix  + ".log2.txt"), "a") as f_hic:
				for k in range(m_log2.shape[0]):
					for l in range(m_log2.shape[1])[k:]:
						x = k * binSize
						y = l * binSize
						value = m_log2[k,l]
						if value == value :
							f_hic.write(chrom+"\t"+str(x)+"\t"+str(x+binSize) + "\t" + chrom+"\t"+str(y) +"\t"+str(y+binSize)+ "\t" + str(value)  + "\n" )
