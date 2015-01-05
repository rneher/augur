# create phylogeny from alignment
# raxml will complain about identical sequences, but still includes them all in the resulting tree
# writes out newick.tree file

import os, re, time, glob, shutil
import subprocess
import dendropy
from io_util import *
							
OUTGROUP = 'A/Beijing/32/1992'							
RAXML_LIMIT = 1.0 # in hours			
raxml_bin = '/ebio/ag-neher/home/rneher/software/standard-RAxML/raxmlHPC-SSE3'
							
def cleanup():
	for file in glob.glob("RAxML_*"):
		try:
			os.remove(file)
		except OSError:
			pass    	
		
def delimit_newick(infile_name, outfile_name):
	with open(infile_name, 'r') as file:
		newick = file.read().replace('\n', '')	
		newick = re.sub(r'(A/[^\:^,]+)', r"'\1'", newick)
	with open(outfile_name, 'w') as file:
		file.write(newick)	
									
def main():

	print "--- Tree infer at " + time.strftime("%H:%M:%S") + " ---"

	cleanup()
	viruses = read_json('data/virus_clean.json')
	write_fasta(viruses, 'temp.fasta')

	print "Building initial tree with FastTree"
	os.system("fasttree_imac -gtr -nt -gamma -nosupport -mlacc 2 -slownni temp.fasta > initial_tree.newick")
	delimit_newick("initial_tree.newick", "temp.newick")
	tree = dendropy.Tree.get_from_path("temp.newick", "newick")	
	tree.resolve_polytomies()
	tree.write_to_path("initial_tree.newick", "newick")

	print "RAxML tree optimization with time limit " + str(RAXML_LIMIT) + " hours"
	os.system("/usr/local/EPD/bin/seqmagick convert temp.fasta temp.phyx")
	# using exec to be able to kill process
	end_time = time.time() + int(RAXML_LIMIT*3600)		
	process = subprocess.Popen("exec ",raxml_bin," -f d -T 6 -j -s temp.phyx -n topology -c 25 -m GTRCAT -p 344312987 -t initial_tree.newick", shell=True)
	while (time.time() < end_time):
		if os.path.isfile('raxml_result.topology'):
			break
		time.sleep(30)
	process.terminate()

	checkpoint_files = [file for file in glob.glob("RAxML_checkpoint*")]
	if os.path.isfile('raxml_result.topology'):
		checkpoint_files.append('raxml_result.topology')
	if len(checkpoint_files) > 0:
		last_tree_file = checkpoint_files[-1]	
		shutil.copy(last_tree_file, 'raxml_tree.newick')
	else:
		shutil.copy("initial_tree.newick", 'raxml_tree.newick')

	print "RAxML branch length optimization and rooting"
	os.system(raxml_bin+" -f e -T 6 -s temp.phyx -n branches -c 25 -m GTRGAMMA -p 344312987 -t raxml_tree.newick -o " + OUTGROUP)
	os.rename('RAxML_result.branches', 'data/tree_branches.newick')

	print "RAxML ancestral state inference"
	os.system(raxml_bin+" -f A -T 6 -s temp.phyx -n states -c 25 -m GTRGAMMA -p 344312987 -t data/tree_branches.newick")
	os.rename('RAxML_nodeLabelledRootedTree.states', 'data/tree_states.newick')
	os.rename('RAxML_marginalAncestralStates.states', 'data/states.txt')
	cleanup()	

if __name__ == "__main__":
    main()
