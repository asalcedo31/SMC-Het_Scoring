import numpy as np
import pprint
import pandas as pd
from ccm_ad_flexible import *
import pdb

def tree_id(tree_df):
#generates a single identifier for each tree structure
	i = 0
	sum = 0
	for n in np.flip(tree_df[:,1],0):
		sum += n*(10**i)
		i += 1
	return sum

def linear_2clust():
	tree = Tree(name=1, size=600)
	node2 = tree.create_node(2,tree.root,600)
	return tree

def linear_3clust():
	tree = Tree(name=1, size=400)
	node2 = tree.create_node(2,tree.root,400)
	node3 = tree.create_node(3,node2,400)	
	# node4 = tree.create_node(4,node3,3)	
	# node5 = tree.create_node(5,node4,3)	
	# node6 = tree.create_node(6,node5,3)	
	return tree

def linear_4clust():
	tree = Tree(name=1, size=300)
	node2 = tree.create_node(2,tree.root,300)
	node3 = tree.create_node(3,node2,300)	
	node4 = tree.create_node(4,node3,300)	
	# node5 = tree.create_node(5,node4,3)	
	# node6 = tree.create_node(6,node5,3)	
	return tree

def linear_5clust():
	tree = Tree(name=1, size=240)
	node2 = tree.create_node(2,tree.root,240)
	node3 = tree.create_node(3,node2,240)	
	node4 = tree.create_node(4,node3,240)	
	node5 = tree.create_node(5,node4,240)	
	# node6 = tree.create_node(6,node5,3)	
	return tree

def linear_6clust():
	tree = Tree(name=1, size=200)
	node2 = tree.create_node(2,tree.root,200)
	node3 = tree.create_node(3,node2,200)	
	node4 = tree.create_node(4,node3,200)	
	node5 = tree.create_node(5,node4,200)	
	node6 = tree.create_node(6,node5,200)	
	return tree

def linear_3clust_uneven(prop=[0.25,0.25,0.5],nssm=1200):
	nssm = float(nssm)
	tree = Tree(name=1, size=nssm*prop[0])
	node2 = tree.create_node(2,tree.root,size=nssm*prop[1])
	node3 = tree.create_node(3,node2,size=nssm*prop[2])	
	return tree

def linear_4clust_uneven(prop=[0.2,0.4, 0.2, 0.2],nssm=1200):
	nssm = float(nssm)
	tree = Tree(name=1, size=nssm*prop[0])
	node2 = tree.create_node(2,tree.root,size=nssm*prop[1])
	node3 = tree.create_node(3,node2,size=nssm*prop[2])	
	node4 = tree.create_node(4,node3,size=nssm*prop[3])	
	return tree

def linear_5clust_uneven(prop=[0.2,0.1, 0.2, 0.2,0.1],nssm=4000):
	nssm = float(nssm)
	tree = Tree(name=1, size=nssm*prop[0])
	node2 = tree.create_node(2,tree.root,size=nssm*prop[1])
	node3 = tree.create_node(3,node2,size=nssm*prop[2])	
	node4 = tree.create_node(4,node3,size=nssm*prop[3])	
	node5 = tree.create_node(5,node4,size=nssm*prop[4])	
	return tree



def permute_tree(tree, tree_dict):
	print tree.tree_struct()
	#starting at the bottom most node switch parents to other nodes with a lower number
	#baseline tree
	tree_array = tree.tree_struct()
	for i in range(len(tree.nodes),1,-1):
		#start with the bottom most node
		node = tree.get_node(i)
	
		#recurse through the nodes above it
		for j in range(1,i):

			#quickly check if the tree you get after switching nodes has been seen before
			new_parent = tree.get_node(j)
			new_tree = tree_array
			new_tree[i-1,1] = j
			new_tree_id = tree_id(new_tree)
			#if the proposed parent only has at most one child and the tree is new then make ti
			if( len(new_parent.children) < 2 and new_tree_id not in tree_dict):
				tree.switch_parent(i,j)
				tree.standard_node_naming()
				tree_df = tree.tree_struct(plot_labels=True,uniform=True)

				#give the tree a unique identifier and store in dict
				tree_num = tree_id(tree_df)
				
				# if the tree is new add it to the dict and recurse
				if tree_num not in tree_dict:
					tree_dict[tree_num] = tree_df[:,1]
					permute_tree(tree,tree_dict)

def output_trees(tree_dict=None, print_trees=False, prefix=None, sizes = None):
	clusters = np.arange(len(tree_dict.values()[0]))+1
	print "clusters", clusters
	i = 0
	if(prefix is None and print_trees == True):
		prefix = 'tree_'+ str(len(clusters)) + "clust"
			
	for parents in tree_dict.values():
		tree_df = np.vstack((clusters,parents))
		tree_df = np.transpose(tree_df)
		tree = tree_from_df(tree_df, sizes=sizes)
		tree_mistakes(tree, prefix + "_"+str(i)+"_")	
		if (print_trees == True):
			outfile = prefix + "_" + str(i) + ".txt"
			print outfile
			print i
			tree.out_1C(prefix + "_" + str(i) + "_1C.txt")
			np.savetxt(outfile,tree_df, fmt='%1i' ,delimiter="\t")
		i += 1

	return (tree_df)

def new_tree_mistakes(tree_dict, i,sizes=None):
	clusters = np.arange(len(tree_dict.values()[0]))+1
	# for parents in tree_dict.values()[0]:
	parents = tree_dict.values()[i]
	
	print(sizes)
	tree_df = np.vstack((clusters,parents))
	tree_df = np.transpose(tree_df)
	# print tree_df
	tree = tree_from_df(tree_df, sizes=sizes)
	tree.tree_struct(plot_labels=True)
	tree.out_1C()
	linearize(tree, "unique_6clust_"+str(i)+"_")


if __name__ == '__main__':


###Standard trees###

	# trees = [ linear_2clust, linear_3clust, linear_4clust, linear_5clust, linear_6clust]

	# i = 2
	# for tree_func in trees:
	# 	tree = tree_func()
	# 	tree_dict = dict()
	# 	sizes = tree.out_1C()['sizes']
	# 	permute_tree(tree,tree_dict)
	# 	pprint.pprint(tree_dict)
	# 	output_trees(tree_dict, print_trees= True, prefix="unique_" + str(i) + "clust", sizes=sizes)
	# 	i += 1

##titration trees###
	proportions = [[0.25,0.25,0.5],
						[0.25,0.5,0.25],
						[0.5,0.25,0.25],
						]
	for i in range(3):
		tree = linear_3clust_uneven(proportions[i])
		tree.tree_struct(outfile="unique_3clust_ssm_"+str(i)+ ".txt",plot_labels=True)
		tree.out_1C("unique_3clust_ssm_" +str(i)+"_1C.txt")
		tree_mistakes(tree, "unique_3clust_ssm_"+ str(i)+ "_")

	proportions = [[0.2,0.2,0.2,0.4],
						[0.2,0.4,0.2,0.2],
						[0.4,0.2,0.2,0.2]						
						]
	for i in range(3):
		tree = linear_4clust_uneven(proportions[i])
		tree.tree_struct(outfile="unique_4clust_ssm_"+str(i)+ ".txt",plot_labels=True)
		tree.out_1C("unique_4clust_ssm_" +str(i)+"_1C.txt")
		tree_mistakes(tree, "unique_4clust_ssm_"+ str(i)+ "_")


###monoclonal tree

	tree = Tree(name=1, size=1200)
	tree.tree_struct(outfile="monoclonal_truth_3A.txt",plot_labels=True)
	tree.out_1C(outfile="monoclonal_truth_1C.txt")
	tree.ccm(outfile="monoclonal_truth_2A.txt")
	split_bottom(tree,"monoclonal_")
	add_extra_bottom(tree,"monoclonal_")

###debug one tree and one scenario###
	# print "linear truth"
	# test = np.asarray([[1,2,3],[0,1,2]])
	# tree = tree_from_df(np.transpose(test),sizes=np.repeat(20,3))
	# print tree.ccm(outfile="truth_2A.txt")

	# tree.out_1C()
	# print tree.tree_struct(plot_labels=True, tier=True)
	# merge_top_and_add_extra(tree,"")
	# tree.tree_struct(plot_labels=True)

	# test_all()

	# print "standard"
	# print tree.standard_node_naming(as_plot_names = False)
	# print tree.tree_struct(plot_labels=True, tier=True)
	# tree.assign_tiers()
	# node4 = tree.get_node(3)
	# tree.gather_within_range(node4,5)
	# orig_1C = tree.out_1C()
	# tree.extra_node(0.75, 3,new_node_name=9, max_dist=3,all_nodes=False)
	# add_intermediate_extra_bottom(tree,"")
	# merge_top_and_add_extra(tree,"")
	# merge_top_and_add_extra(tree,"")
	# linearize(tree,"")
	# tree.tree_struct(plot_labels=True)
	# new_1C = tree.out_1C()
	# print sum(orig_1C['sizes']) == sum(new_1C['sizes'])
	# print sum(orig_1C['sizes'])
	# split_tree = split_bottom(tree,"")
	# tree.ccm()
	# split_tree.ccm(output2A=True, outfile="test2A.txt")

	# sizes = tree.out_1C()['sizes']
	# print(sizes)
	# tree_dict = dict()
	# tree.assign_tiers()
	# permute_tree(tree,tree_dict)
	# pprint.pprint(tree_dict)
	# new_tree_mistakes(tree_dict,1,sizes)

	# tree_df = tree.tree_struct(tier = True, plot_labels=True)

	# tree = linear_4clust()
	# orig_tree_df = tree.tree_struct()
	# orig_tree_num = tree_id(orig_tree_df)
	# print "orig"
	# print tree.tree_struct()
	# tree_dict = dict()
	# tree_dict[orig_tree_num] = orig_tree_df[:,1]
	# permute_tree(tree,tree_dict)
	# pprint.pprint(tree_dict)
	# output_trees(tree_dict, print_trees= True, prefix="unique_2clust")


	
	# tree.out_1C("test_1C.txt")
	# tree_dict = dict()
	# permute_tree(tree,tree_dict)
	# output_trees(tree, print_trees= True, prefix="unique_3clust_ssm")


	# for i in range(6,7):
	# 	new_tree_mistakes(tree_dict,i)



	# tree.uniform_node_naming()
	# tree_df = tree.tree_struct(plot_labels = True,uniform=True)
	# print tree_id(tree_df)

	# tree = branching_test4()
	# tree_df = tree.tree_struct(plot_labels = True)
	# print tree_df

	# tree.uniform_node_naming()
	# tree_df = tree.tree_struct(plot_labels = True,uniform=True)
	# print tree_id(tree_df)



