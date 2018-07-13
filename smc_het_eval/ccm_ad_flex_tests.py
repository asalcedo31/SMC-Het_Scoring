
def linear_truth():
	#truth1
	tree = Tree(name='N1', size=3)
	node2 = tree.create_node('N2',tree.root,2)
	node3 = tree.create_node('N3',node2,3)
	node4 = tree.create_node('N4',node3,3)
	return tree

def linear_branching_truth():
	#truth2
	tree = Tree(name='N1', size=3)
	node2 = tree.create_node('N2',tree.root,2)
	node3 = tree.create_node('N3',node2,3)
	node4 = tree.create_node('N4',node2,3)
	return tree



def branching_to_linear(tree):
	#truth1 as a mistake for truth1
	tree.switch_parent('N4','N3')
	return tree

def split_four(tree, prop_split=0.5):
	#corresponds to 2
	# print tree.nodes[1].children
	if(tree.nodes[1].children.size > 1): #if branching
		tree.switch_parent('N4','N3')
	tree.split_node('N4','N5',prop_split=prop_split, same=True)
	return tree

def collapse_bottom_clusters(tree):
	#corresponds to 3
	tree.merge_nodes('N3','N4')
	return tree

def collapse_all_clusters_and_extra(tree, extra_prop=0.5):
	#corresponds to 15
	tree.collapse_node('N2')
	tree.extra_node(extra_prop, 'N2')
	return tree

def extra_two(tree, extra_prop=0.5):
	#corresponds to 4
	tree.extra_node(extra_prop,'N2')
	return tree

def collapse_mid(tree):
	#corresponds to 5
	tree.merge_nodes('N2','N3')
	return(tree)


def extra_four(tree, extra_prop=0.5):
	#corresponds to 6
	tree.extra_node(extra_prop,'N4')
	return tree

def split_four_lin(tree):
	#corresponds to 13
	if(tree.nodes[1].children.size > 1): #if branching
		tree.switch_parent('N4','N3')
	tree.split_node('N4','N5',same=False)
	return tree

def split_three_from_two(tree):
	#corresponds to 13
	if(tree.nodes[1].children.size > 1): #if branching
		tree.switch_parent('N4','N3')
	tree.split_node('N3','N5',same=True)
	return tree

def collapse_all_bottom(tree):
	tree.collapse_node('N2')
	return tree


def extra_one(tree, extra_prop=0.5):
	#corresponds to 11
	tree.extra_node(extra_prop,'N1', transfer_children=True)
	return tree

def extra_root(tree, extra_prop=0.5):
	#corresponds to 12
	tree.extra_node(extra_prop,None, transfer_children=True)
	return tree

def extra_one_switch_three(tree, extra_prop=0.5):
	#corresponds to 8
	if(tree.nodes[1].children.size > 1): #if branching
		tree.switch_parent('N4','N3')	
	tree.switch_parent('N3', 'N1')
	tree.extra_node(extra_prop,'N1', transfer_children=True)
	return tree

def three_from_one(tree):
#corresponds to 9
	if(tree.nodes[1].children.size > 1): #if branching
		tree.switch_parent('N4','N3')	
	tree.switch_parent('N3','N1')
	return(tree)

def three_from_one(tree):
#corresponds to 10
	if(tree.nodes[1].children.size == 1): #if branching
		tree.switch_parent('N4','N2')	
	tree.switch_parent('N3','N1')
	return(tree)


scenarios ={
	'T1' : linear_truth,
	'T2' : linear_branching_truth,
	'S2' : split_four,
	'S3' : collapse_bottom_clusters,
	'S4' : extra_two
}
