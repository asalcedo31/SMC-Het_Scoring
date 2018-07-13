import numpy as np
import pprint
import pandas as pd
from ccm_ad_flexible import *
import pdb

def tree_id(tree_df):
	i = 0
	sum = 0
	for n in np.flip(tree_df[:,1],0):
		# print "n", n, "i", i,n*(10**i)
		sum += n*(10**i)
		i += 1
	return sum

def linear_4clust():
	tree = Tree(name=1, size=3)
	node2 = tree.create_node(2,tree.root,2)
	node3 = tree.create_node(3,node2,3)	
	node4 = tree.create_node(4,node3,3)	
	node5 = tree.create_node(5,node4,3)	
	node6 = tree.create_node(6,node5,3)	
	return tree



def permute_tree(tree, tree_dict):
	print tree.tree_struct()
	#starting at the bottom most node switch parents to other nodes with a lower number
	#baseline tree
	tree_array = tree.tree_struct()
	for i in range(len(tree.nodes),1,-1):
		#start with the bottom most node
		node = tree.get_node(i)
		# print "i", i, "node", node.name
	
		#recurse through the nodes above it
		for j in range(1,i):

			# print "j",j
			#quickly check if the tree you get after switching nodes has been seen before
			new_parent = tree.get_node(j)
			new_tree = tree_array
			new_tree[i-1,1] = j
			new_tree_id = tree_id(new_tree)
			#if the proposed parent only has at most one child and the tree is new then make ti
			
			
			if( len(new_parent.children) < 2 and new_tree_id not in tree_dict):
				print "orig"
				print tree.tree_struct(plot_labels=True)
	

				tree.switch_parent(i,j)
				# print"orig"
				print tree.tree_struct()
				tree.standard_node_naming()
				tree_df = tree.tree_struct(plot_labels=True,uniform=True)

				#give the tree a unique identifier and store in dict
				tree_num = tree_id(tree_df)
				

				# if the tree is new add it to the dict and recurse
				if tree_num not in tree_dict:
					# if ( tree_num==11334):
					# 	pdb.set_trace()
					# 	tree.standard_node_naming()
				
					print "new tree ", i, j 
					print tree.tree_struct(plot_labels=True, uniform = True, tier=True)				
					print "tree id", tree_num, "tiers", tree.tiers
					tree_dict[tree_num] = tree_df[:,1]
					permute_tree(tree,tree_dict)

def output_trees(tree_dict=None, print_trees=False, prefix=None):
	clusters = np.arange(len(tree_dict.values()[0]))+1
	print "clusters", clusters
	i = 0
	if(prefix is None and print_trees == True):
		prefix = 'tree_'+ str(len(clusters)) + "clust"
			
	for parents in tree_dict.values():
		tree_df = np.vstack((clusters,parents))
		tree_df = np.transpose(tree_df)

		if (print_trees == True):
			outfile = prefix + "_" + str(i) + ".txt"
			print outfile
			print i
			np.savetxt(outfile,tree_df, fmt='%1i')
		i += 1

	return (tree_df)

def new_tree_mistakes(tree_dict, i):
	clusters = np.arange(len(tree_dict.values()[0]))+1
	# for parents in tree_dict.values()[0]:
	parents = tree_dict.values()[i]
	
	tree_df = np.vstack((clusters,parents))
	tree_df = np.transpose(tree_df)
	# print tree_df
	tree = tree_from_df(tree_df)
	tree_mistakes(tree)


if __name__ == '__main__':

	tree = linear_4clust()
	orig_tree_df = tree.tree_struct()
	orig_tree_num = tree_id(orig_tree_df)
	print "orig"
	print tree.tree_struct()
	tree_dict = dict()
	tree_dict[orig_tree_num] = orig_tree_df[:,1]
	permute_tree(tree,tree_dict)
	pprint.pprint(tree_dict)
	# output_trees(tree_dict, print_trees= True, prefix="unique_6clust")
	for i in range(3):
		new_tree_mistakes(tree_dict,i)

	# print "branching truth"
	# tree = branching_test5()
	# tree.assign_tiers()
	# tree_df = tree.tree_struct(tier = True)
	# print tree_df
	# print tree.tiers

	# tree.uniform_node_naming()
	# tree_df = tree.tree_struct(plot_labels = True,uniform=True)
	# print tree_id(tree_df)

	# tree = branching_test4()
	# tree_df = tree.tree_struct(plot_labels = True)
	# print tree_df

	# tree.uniform_node_naming()
	# tree_df = tree.tree_struct(plot_labels = True,uniform=True)
	# print tree_id(tree_df)



