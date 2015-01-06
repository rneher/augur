# clean, reroot, ladderize newick tree
# output to tree.json

import os, re, time, datetime
import dendropy
from io_util import *
from tree_LBI import *

OUTGROUP = 'A/Beijing/32/1992'
# cluster position in numbering starting at 0
cluster_positions = np.array(sorted([188,192,155,158,157, 154, 144]), dtype=int)+16


def delimit_newick(infile_name):
	with open(infile_name, 'r') as file:
		newick = file.read().replace('\n', '')	
		newick = re.sub(r'(A/[^\:^,^)]+)', r"'\1'", newick)
	return newick
		
def crossref_import(branches_tree_file, states_tree_file, states_file):
	"""RAxML won't return a single NEWICK tree with both ancestral states and branch lengths"""
	"""This opens the necessary RAxL output files and outputs a single Dendropy tree"""
	label_to_seq = {}
	with open(states_file) as file:
		for line in file:
			(label, seq) = line.split()
			label_to_seq[label] = seq
	branches_tree = dendropy.Tree.get_from_string(delimit_newick(branches_tree_file), "newick", as_rooted=True)
	states_tree = dendropy.Tree.get_from_string(delimit_newick(states_tree_file), "newick", as_rooted=True)	
	for (bn, sn) in zip(branches_tree.postorder_node_iter(), states_tree.postorder_node_iter()):
		if sn.label:
			bn.seq = label_to_seq[sn.label]
	return branches_tree

def date_to_str(date):
	if isinstance(date, datetime.date):
		return date.isoformat()
	elif isinstance(date,basestring):
		return date
	elif isinstance(date, int) or isinstance(date, float):
		return datetime.date.fromordinal(int(date)).isoformat()
	else:
		print "unknown date format", date
		return ""

def to_json(node):
	json = {}
	for attr, val in node.__dict__.iteritems():
		if type(val) in [float, int, np.float64]:
			json[attr] = round(val,5)
		elif isinstance(val, basestring):
			json[attr] = val.replace("'",'')

	if hasattr(node, 'clade'):
		json['clade'] = node.clade
	if node.taxon:
		json['strain'] = str(node.taxon).replace("'", '')
	if hasattr(node, 'yvalue'):
		json['yvalue'] = round(node.yvalue, 5)
	if hasattr(node, 'xvalue'):
		json['xvalue'] = round(node.xvalue, 5)
	if hasattr(node, 'date'):
		json['date'] = date_to_str(node.date)
	if hasattr(node, 'seq'):
		json['seq'] = node.seq
	if hasattr(node, 'koel'):
		json['koel'] = node.koel
	if node.child_nodes():
		json["children"] = []
	for ch in node.child_nodes():
		json["children"].append(to_json(ch))
	return json

#def to_json(node):
#	import numpy as np
#	json = {}
#	for attr, val in node.__dict__.iteritems():
#		if type(val) in [float, int, np.float64]:
#			json[attr] = round(val,5)
#		elif isinstance(val, basestring):
#			json[attr] = val.replace("'",'')
#	if node.taxon:
#		json['strain'] = str(node.taxon).replace("'", '')
#	if node.child_nodes():
#		json["children"] = []
#		for ch in node.child_nodes():
#			json["children"].append(to_json(ch))
#	return json
	
def get_yvalue(node):
	"""Return y location based on recursive mean of daughter locations"""	
	if hasattr(node, 'yvalue'):
		return node.yvalue	
	if node.child_nodes():
		mean = 0
		for ch in node.child_nodes():
			mean += get_yvalue(ch)
		return mean / float(len(node.child_nodes()))
		
def get_xvalue(node):
	"""Return x location based on total distance from root"""
	root = node.get_tree_root()
	return node.get_distance(root)

def remove_outgroup(tree):
	"""Reroot tree to outgroup"""
	outgroup_node = tree.find_node_with_taxon_label(OUTGROUP)	
	if outgroup_node:
#		tree.to_outgroup_position(outgroup_node, update_splits=False)
		tree.prune_subtree(outgroup_node)
	else:
		print "outgroup not found"
def collapse(tree):
	"""Collapse short edges to polytomies"""
	for edge in tree.postorder_edge_iter():
		if edge.length < 0.00001 and edge.is_internal():
			edge.collapse()
			
def reduce(tree):
	"""Remove outlier tips"""
	for node in tree.postorder_node_iter():
		if node.edge_length > 0.04 and node.is_leaf():
			parent = node.parent_node
			parent.remove_child(node)
			
def ladderize(tree):
	"""Sorts child nodes in terms of the length of subtending branches each child node has"""
	node_desc_counts = {}
	for node in tree.postorder_node_iter():
		if len(node._child_nodes) == 0:
			node_desc_counts[node] = node.edge_length
		else:
			total = 0
			if node.edge_length > 0:
				total += node.edge_length			
			for child in node._child_nodes:
				total += node_desc_counts[child]
			node_desc_counts[node] = total
			node._child_nodes.sort(key=lambda n: node_desc_counts[n], reverse=True)			

def add_node_attributes(tree):
	"""Add clade, xvalue and yvalue attributes to all nodes in tree"""
	clade = 0
	yvalue = 0
	for node in tree.postorder_node_iter():
		node.clade = clade
		clade += 1
		if node.is_leaf():
			node.yvalue = yvalue
			yvalue += 1

	for node in tree.postorder_node_iter():
		node.yvalue = get_yvalue(node)
		node.xvalue = node.distance_from_root()
		
def layout(tree):
	"""Set yvalue of tips by post-order traversal"""
	yvalue = 0	
	distance_matrix = dendropy.treecalc.PatristicDistanceMatrix(tree)
	tips = [node for node in tree.leaf_iter()]
	tips[0].yvalue = yvalue
	for (a,b) in zip(tips[:-1], tips[1:]):
		d = distance_matrix(a.taxon, b.taxon)
	#	print str(a.taxon) + " to " + str(b.taxon) + ": " + str(d)
		if b.is_leaf():
			yvalue += d
			b.yvalue = yvalue
			
	for node in tree.postorder_node_iter():
		node.yvalue = get_yvalue(node)		

def add_virus_attributes(viruses, tree):
	"""Add date and seq attributes to all tips in tree"""
	strain_to_date = {}
	strain_to_seq = {}
	for v in viruses:
		strain_to_date[v['strain']] = v['date']
		strain_to_seq[v['strain']] = v['seq']
	for node in tree.postorder_node_iter():
		strain = str(node.taxon).replace("'", '')
		if strain_to_date.has_key(strain):
			node.date = strain_to_date[strain]
		if strain_to_seq.has_key(strain):
			node.seq = strain_to_seq[strain]


def add_LBI(tree):
	print "calculate local branching index"
	T2 = get_average_T2(tree, 365)
	tau =  T2*2**-4
	print "avg pairwise distance:", T2
	print "memory time scale:", tau
	calc_LBI(tree, tau = tau)

def add_Koel_gt(tree):
	for node in tree.postorder_node_iter():
		node.koel = Koel_gt(node.seq);

def coordinates_to_region(lat, lng):
    '''
    returns the region based on the geographic sectors defined by longitude and latitude
    argument:
    lat  -- latitude
    lng  -- longitude
    '''
    if lat>0:
        if lng<-20:
            return 'north_america'
        elif lng>50:
            return 'asia'
        else:
            return 'europe'
    else:
        if lng<-20:
            return 'south_america'
        else:
            return 'oceania'

def add_place(places_to_coordinates, place):
    from geopy import geocoders
    g = geocoders.GoogleV3()
    g.timeout=10
    loc = g.geocode(place.replace('_', ' '))
    time.sleep(0.2)
    print place, loc
    try:
        if loc:
            country = loc[0].split(',')[-1].strip()
            country = country.encode('ascii', 'replace')
            location = loc[1]
            places_to_coordinates[place]={}
            places_to_coordinates[place]['country']=country
            places_to_coordinates[place]['lat']=loc[0]
            places_to_coordinates[place]['lng']=loc[1]
            print place, loc
        else:
            print "place not resolved"
    except:
        print "ERROR",place

def add_geo_info(tree):
    import pickle
    with open("data/flubase_places.pickle","r") as pfile:
        places_to_coordinates = pickle.load(pfile)
    for node in tree.leaf_iter():
        try:
            place = node.strain.split('/')[1]
        except:
            place = node.taxon.label.split('/')[0]
        if place not in places_to_coordinates:
            add_place(places_to_coordinates, place)
        if place in places_to_coordinates:
            node.country = places_to_coordinates[place]['country'].upper()
            node.lat =places_to_coordinates[place]['lat']
            node.lng =places_to_coordinates[place]['lng']
            node.region = coordinates_to_region(node.lat, node.lng)
        else:
            print "no geo info for",place
            node.country = "undefined"
            node.lat = "undefined"
            node.lng = "undefined"
            node.region = "undefined"

    with open("data/flubase_places.pickle","w") as pfile:
        pickle.dump(places_to_coordinates,pfile)

def Koel_gt(seq):
	import Bio
	aaseq = str(Bio.Seq.Seq(seq).translate())[:-1]
	if '*' not in aaseq:
		return "".join([aaseq[i] for i in cluster_positions])
	else:
		print "translation has stop"
		return ''

def main():

	print "--- Tree clean at " + time.strftime("%H:%M:%S") + " ---"
		
	viruses = read_json('data/virus_clean.json')
	tree = crossref_import('data/tree_branches.newick', 'data/tree_states.newick', 'data/states.txt')
	print "Remove outgroup"
	remove_outgroup(tree)
	print "Remove outlier branches"	
	reduce(tree)
	print "Collapse internal nodes"		
	collapse(tree)	
	print "Ladderize tree"	
	ladderize(tree)
	print "added attributes"
	add_node_attributes(tree)
	add_virus_attributes(viruses, tree)
	print "add Koel genotypes"
	add_Koel_gt(tree)
	print "add geo info"
	add_geo_info(tree)
	#add_LBI(tree)
	write_json(to_json(tree.seed_node.child_nodes()[0]), "data/tree_clean.json")
	return tree

if __name__ == "__main__":
	tree = main()
