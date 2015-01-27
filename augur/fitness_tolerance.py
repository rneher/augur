import numpy as np
from Bio import Seq, AlignIO
import time
from io_util import read_json
from io_util import write_json

def load_mutational_tolerance():
	fname = '../data/Thyagarajan_Bloom_HA_fitness.txt'
	with open(fname) as f:
		aa = map(lambda x:x.split('_')[1], f.readline().strip().split()[3:])
	sites = np.loadtxt(fname, usecols=[0], dtype=int)
	wt_aa = np.loadtxt(fname, usecols=[1], dtype='S1')
	aa_prob = np.loadtxt(fname, usecols=range(3,23), dtype=float)
	return aa, sites, wt_aa, aa_prob

def calc_fitness_tolerance(aa_seq, aa_prob, aa, indices, beta = 1.0):
	'''
	determine the indices of aligned amino acids and sum the logged probabilities
	'''
	H3_aa_indices = [aa.index(aa_seq[p]) if aa_seq[p] in aa else -1 for p in indices]
	return np.sum(np.log(aa_prob[(np.arange(len(indices)), H3_aa_indices)])*beta)

def assign_fitness(viruses):
	'''
	loops over all viruses, translates their sequences and calculates the virus fitness
	'''
	aa, sites, wt_aa, aa_prob = load_mutational_tolerance()
	aln = AlignIO.read('../data/H1_H3.fasta', 'fasta')
	# returns true whenever either of the sequences have a gap
	aligned = (np.array(aln)!='-').min(axis=0)
	# map alignment positions to sequence positions, subset to aligned amino acids
	indices = {}
	for seq in aln:
		indices[seq.name] = (np.cumsum(np.fromstring(str(seq.seq), dtype='S1')!='-')-1)[aligned]

	# make a reduced set of amino-acid probabilities that only contains aligned positions
	aa_prob=aa_prob[indices['H1'],:]
	# attach another column for non-canonical amino acids
	aa_prob = np.hstack((aa_prob, 1e-5*np.ones((aa_prob.shape[0],1))))
	for virus in viruses:
		virus['fitness_tolerance'] = calc_fitness_tolerance(Seq.translate(virus['seq']), 
															aa_prob, aa, indices['H3'])

def main():

	print "--- Mutational tolerance at " + time.strftime("%H:%M:%S") + " ---"
	viruses = read_json('data/virus_clean.json')
	assign_fitness(viruses)
	write_json(viruses, "data/virus_tolerance.json")

if __name__=='__main__':
	main()
