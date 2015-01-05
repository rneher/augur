import time, os
import virus_filter, virus_align, virus_clean
import tree_infer, tree_clean, tree_frequency, tree_auspice

def main(y=3, vpm=75):
	"""Run full pipeline"""
	
	print "--- Start processing at " + time.strftime("%H:%M:%S") + " ---"	
    
	virus_filter.main(y,vpm)			# Filter sequences
	virus_align.main()			# Align sequences
	virus_clean.main()			# Clean sequences	
	tree_infer.main()			# Make tree
	tree_clean.main()			# Clean tree	
#	tree_frequency.main()		# Add clade frequencies
	tree_auspice.main()			# Streamline tree for auspice

if __name__ == "__main__":
    main()
