import numpy as np;
from make_ad_nvar import *
from make_ccm import *
import pprint
import pandas as pd
import operator
import pdb
import copy
import re
from ccm_ad_flex_tests import *


class Node:
	def __init__(self,value, name, parent=None):
		self.value = value
		self.children = np.array([])
		self.desc = np.array([])
		self.parent = parent
		self.name = name
		self.plot_name = name
		self.desc_depth = 0
		self.tier = 0
		self.desc_idx = np.array([],np.int32)
		if parent is not None:
			parent.add_children([self])
		else:
			self.tier = 1

	def add_children(self, children):
		#children are node instances
		#put single nodes into a list
		if isinstance(children, list)==False and isinstance(children,np.ndarray) == False:
			children = [children]
		if self.children is not None:
			for child in children:
				self.children = np.append(self.children,child)
	def remove_child(self, child_names):
		#child_names is a string or string list
		#put single node names into a list 
		if isinstance(child_names, basestring) or isinstance(child_names, int):
			child_names = [child_names]
		#list of indices to remove
		to_delete = np.array([],np.int32)
		for idx in range(len(self.children)):
			if self.children[idx].name in child_names:
				to_delete = np.append(to_delete, idx)
		self.children = np.delete(self.children, to_delete)
		return None
	def gather_desc_nodes(self, top_node= None):
		if top_node is None:
			top_node = self
			top_node.desc = np.array([])
		if (len(self.children) > 0 ):
			top_node.desc_depth += 1
		for child in self.children:
			top_node.desc = np.append(top_node.desc, child.name)
			if len(child.children) > 0:
				#get the names of the descendants of each child
				child.gather_desc_nodes(top_node)

	def gather_child_idx(self,top_node=None):
		#top node is the node for which we will gather all of the descendant indices
		if top_node is None:
			top_node = self
		for child in self.children:
			#get the indices for each child

			top_node.desc_idx = np.append(top_node.desc_idx, child.value)
			if len(child.children) > 0:
				#get the indices for the descendants of each child
				child.gather_child_idx(top_node)
	def collapse_children(self, full=True):
		#ALL of the node's descendants (not just children) are merged into the cluster
		if full == True:
			top_node = self
			top_node.desc_idx = np.array([], np.int32)
			self.gather_child_idx(top_node)
			#all of the node's descendants now become part of the cluster
			self.value = np.append(self.value,self.desc_idx)
			self.desc_idx = np.array([],np.int32)
			self.children = np.array([])

class Tree:
	def __init__(self,name, size=None, values=None):
		self.nodes = []
		if values is None:
			root_idx = self.node_indices(size)
		else:
			root_idx = values
		self.root = Node(root_idx,name)
		self.nodes = [self.root]
		self.node_tiers = None
		self.tiers = 1

	def node_indices(self,size):
		#generates the mutation indices for a new node
		idx = np.arange(size)
		nssm = self.comp_nssms() #how many ssms are already in the tree
		idx = idx + nssm #make an array of size as the value of the node
		return idx

	def create_node(self,name,parent=None,size=None, values=None):
		if size is not None and values is None:
			idx = self.node_indices(size)
		if values is not None:
			idx = values
		node = Node(idx, name, parent)
		self.nodes= np.append(self.nodes,node)
		return node

	def comp_nssms(self):
		nssms = reduce(lambda x,y: x + y.value.size,self.nodes ,0)
		return(nssms)

	def ccm(self):
		nssms = self.comp_nssms()
		out = np.zeros((nssms,len(self.nodes)))
		i = 0
		for node in self.nodes:
			out[node.value,i] = 1
			i += 1
		#print out
		ccm = np.dot(out,out.T)
		return ccm

	def ad(self):
		nssms = self.comp_nssms()
		out = np.zeros((nssms,nssms))
		for node in self.nodes:
			#clear descendent idx array first
			node.desc_idx = np.array([],np.int32)
			#repopulate
			node.gather_child_idx()
			desc_idx = node.desc_idx
			rows = []
			if isinstance(node.value, np.int64):
				rows = [node.value]
			else:
				rows = node.value
			out[np.ix_(rows, desc_idx)] = 1
		return out
		
	def tree_struct(self,plot_labels=False, uniform = False, outfile = None, tier = False):
		#creates the 3A output with an optional plot labels cloumn that will be included in the printed outfile
		#TODO add rename nodes to numbers function
		node_names = list()
		parent_names = list()
		node_numbers = list()
		parent_numbers = list()
		node_tiers = list()
		for node in self.nodes:
			node_name = node.name
			if node.parent is not None:
				parent_name = node.parent.name
				parent_names.append(node.parent.plot_name)
			else:
				parent_name = 0
				parent_names.append(0)

			
			if isinstance(node_name, basestring):
				node_match = re.match(r'N\d',node_name)
				if node_match is not None:
					node_number = int(node_name.replace('N',''))
				else:
					node_number = node_name
			else:
				node_number = node_name
			if isinstance(parent_name, basestring):
				parent_match = re.match(r'N\d',parent_name)
				if parent_match is not None:
					parent_number = int(parent_name.replace('N',''))
				else:
					parent_number = parent_name
			else:
				parent_number = parent_name

			node_numbers.append(node_number)
			parent_numbers.append(parent_number)
			node_names.append(node.plot_name)
			node_tiers.append(node.tier)




		if plot_labels == False:
			tree_df = pd.DataFrame({'nodes' : node_numbers, 'parent' : parent_numbers, 'tier': node_tiers})

		if plot_labels == True:
			tree_df =  pd.DataFrame({'nodes' : node_numbers, 'parent' : parent_numbers, 'plot_names' : node_names, 'parent_names': parent_names, 'tier': node_tiers})

		#rename any nodes with string names as increasing integers
		if tree_df.dtypes['nodes'] != 'int64':
			tree_df['nodes'] = tree_df['nodes'].astype(str)
			tree_df['parent'] = tree_df['parent'].astype(str)
			# print tree_df.dtypes

			int_nodes = tree_df[tree_df['nodes'].str.contains(r'^[\d]+$')]
			string_nodes = tree_df[tree_df['nodes'].str.contains(r'[a-zA-Z]')]
			# str_idx = string_nodes.index()
			node_names = dict()
			max_int = int(max(int_nodes['nodes']))
			j = 1
			for i in string_nodes.index:
				node = string_nodes['nodes'][i]
				node_names[node] = max_int + j
				tree_df['nodes'][i] = max_int + j
				j += 1
			
			string_parents = tree_df[tree_df['parent'].str.contains(r'[a-zA-Z]')]
			# str_idx = string_nodes.index()
			for i in string_parents.index:
				
				tree_df['parent'][i] = node_names[tree_df['parent'][i]]


		if outfile is not None:
			tree_df.to_csv(outfile, header=None,index=None,sep="\t")
		if (uniform == True and tier == True):
			tree_df = tree_df.sort_values(by='plot_names')
			print tree_df[['nodes','parent','plot_names', 'parent_names', 'tier']]
			return tree_df.as_matrix(columns=['plot_names', 'parent_names'])
		
		if (uniform == True):
			tree_df = tree_df.sort_values(by='plot_names')
			# print tree_df[['nodes','parent','plot_names', 'parent_names', 'tier']]
			return tree_df.as_matrix(columns=['plot_names', 'parent_names'])
		if (tier == True):
			print tree_df[['nodes','parent','tier']]
		
		if plot_labels == True:
			print tree_df[['nodes','parent','plot_names', 'parent_names']]
		
		return tree_df.as_matrix(columns=['nodes','parent'])


	def get_node(self,name, return_idx=False):
		for idx in range(len(self.nodes)):
			node = self.nodes[idx]
			if node.name == name:
				if return_idx == True:
					return node, idx
				else: 
					return node

	def split_node(self,name,new_name, prop_split=0.5,same=True):
		node = self.get_node(name)
		node_idx = node.value
		prop_keep = 1- prop_split
		num_to_keep = int(prop_keep*len(node_idx))
		if num_to_keep == 1:
			idx_keep = [0]
		else:
			idx_keep = range(num_to_keep)
		idx_split = range(num_to_keep, len(node_idx))
		node.value = node_idx[idx_keep]
		if same == True:
			parent = node.parent
		else:
			parent = node
		new_node = self.create_node(new_name,parent=parent,values=node_idx[idx_split])

	def merge_nodes(self, node_name1, node_name2):
		node1 = self.get_node(node_name1)
		node2, idx2 = self.get_node(node_name2, return_idx=True)
		node1.value = np.concatenate((node1.value,node2.value))
		node1.remove_child(node_name2)
		if len(node2.children) > 0 :
			node1.add_children(node2.children)
		self.nodes = np.delete(self.nodes, idx2)

	def collapse_node(self, node_name):
		node = self.get_node(node_name)
		node.desc = np.array([])
		node.gather_desc_nodes()
		node.collapse_children()
		for name in node.desc:
			d_node, d_idx = self.get_node(name,return_idx = True)
			self.nodes = np.delete(self.nodes, d_idx)

	def _gather_tiered_nodes(self):
		self.node_tiers = dict()
		for tier in range(1,self.tiers+1):
			self.node_tiers[tier] = list()

		for node in self.nodes:
			try:
				self.node_tiers[node.tier].append(node)
			except:
				die
				# pdb.set_trace()

	def _tier_assignment(self, node = None):
		if node is None:
			node = self.root
		children = node.children
		
		for child in children:
			child.tier = child.parent.tier +1
			if child.tier > self.tiers:
				self.tiers = child.tier
			self._tier_assignment(node=child)

	def assign_tiers(self):
		self._tier_assignment()
		self._gather_tiered_nodes()

	def get_tier(self, tier):
		if self.node_tiers is None:
			self.assign_tiers()
		
		return(self.node_tiers[tier])

	

	def standard_node_naming(self):
		self.assign_tiers()
		max_node_num = 1
		for tier in range(1,self.tiers+1):
			# print "tier ", tier
			if tier == 1:
				node = self.root
				node.plot_name = 1
				self.max_node_num = 1
				continue
			
			nodes = self.get_tier(tier)
			
			n_children = map(lambda x: len(x.children),nodes)
			map(lambda x: x.gather_desc_nodes(),nodes)
			n_desc = map(lambda x: len(x.desc),nodes)
			names = map(lambda x: x.name,nodes)
			parent_names = map(lambda x: x.parent.plot_name,nodes)

			# print "nchildren"
			# print n_children
			# print "ndesc"
			# print n_desc
			# print "names"
			# print names

			sorted_nodes = [x for _,_,_,_,x in sorted(zip(n_children, n_desc,parent_names, names, nodes), key = operator.itemgetter(0,1,2))]

			for node in sorted_nodes:
				node.plot_name = max_node_num + 1
				max_node_num += 1 

			# map(lambda x: x.name,nodes)


	def uniform_node_naming(self,node=None):
		if node is None:
			node = self.root
			node.plot_name = 1
			self.max_node_num = 1
		print node.name
		# print "self.max_node_num", self.max_node_num
		children = node.children
		if ( len(children) == 1 ):
			children[0].plot_name = self.max_node_num +1
			self.max_node_num = children[0].plot_name
			self.uniform_node_naming(children[0])
		
		elif (len(children) == 2 ):
			n_children = [len(children[0].children), len(children[1].children)]
			children[0].gather_desc_nodes()
			children[1].gather_desc_nodes()
			# print children
			n_desc = [len(children[0].desc),len(children[1].desc)]
			names = [children[0].name, children[1].name]
			# print children[1].desc_depth
			print n_children
			print n_desc
			if (n_children[0] != n_children[1]):
				#if one has more children it will get the higher node number
				children[np.argmin(n_children)].plot_name = self.max_node_num + 1
				children[np.argmax(n_children)].plot_name = self.max_node_num + 2
				self.max_node_num = children[np.argmax(n_children)].plot_name
				children = [x for _,x in sorted(zip(n_desc,children))]
				

			elif (n_desc[0] != n_desc[1]):
				print "most desc ", children[np.argmax(n_desc)].name
				#if both have the name number of children the one with more descendents gets the higher number
				children[np.argmin(n_desc)].plot_name = self.max_node_num + 1
				children[np.argmax(n_desc)].plot_name = self.max_node_num + 2
				self.max_node_num = children[np.argmax(n_desc)].plot_name
				children = [x for _,x in sorted(zip(n_desc,children))]
				

			else: 
				children = [x for _,x in sorted(zip(names,children))]
				children[0].plot_name = self.max_node_num + 1
				children[1].plot_name = self.max_node_num + 2
				self.max_node_num = children[1].plot_name
				children = [x for _,x in sorted(zip(names,children))]


			for child in children:
				# print "child ", child.name
				self.uniform_node_naming(child)
		
		else: #no children
			return



	def switch_parent(self, node_name,new_parent_name):
		node = self.get_node(node_name)
		new_parent_node = self.get_node(new_parent_name)
		old_parent_node = self.get_node(node.parent.name)
		
		if new_parent_node == old_parent_node:
			return
		node.parent = new_parent_node
		
		new_parent_node.add_children(node)
		old_parent_node.remove_child(node.name)
		#update descendants of the new parent node
		if len(new_parent_node.desc) > 0:
			new_parent_node.gather_desc_nodes()
		if len(new_parent_node.desc_idx) > 0:
			new_parent_node.gather_child_idx()
		#update descendants of the old parent node
		if len(old_parent_node.desc) > 0:
			old_parent_node.gather_desc_nodes()
		if len(old_parent_node.desc_idx) > 0:
			old_parent_node.gather_child_idx()

	def extra_node(self,extra_prop, parent_name, new_node_name='Extra', transfer_children = False, all_nodes=True, num_nodes = 2):
		drawn_nodes = []
		if all_nodes == True:
			drawn_nodes = self.nodes
		else:
			drawn_nodes = np.random.choice(self.nodes, num_nodes, replace = False)
		extra_idx = np.array([],np.int32)
		for node in drawn_nodes:
			num_taken = int(round(extra_prop*len(node.value)))
			taken_idx = []
			if num_taken == 1:
				taken_idx = 0
			elif num_taken == 0:
				taken_idx = []
			else:
				taken_idx = range(num_taken-1)
			taken = node.value[taken_idx]
			extra_idx = np.append(extra_idx, taken)
			node.value = np.delete(node.value,taken_idx)
		if parent_name is not None:
			parent_node = self.get_node(parent_name)
			extra = self.create_node(new_node_name, parent=parent_node,values=extra_idx)
			if transfer_children == True:
				# extra.add_children(parent_node.children)
				# child_names = [child.name for child in parent_node.children[0:parent_node.children.size-2]]
				# print parent_node.children.size
				[self.switch_parent(child.name,new_node_name) for child in parent_node.children[0:parent_node.children.size-1]]
				# print [child.name for child in parent_node.children[0:parent_node.children.size-1]]
				# print parent_node.children
				
				# parent_node.remove_child(child_names)
				# parent_node.add_children(extra)
				# parent_node.gather_desc_nodes()
				# parent_node.gather_child_idx()
				# extra.gather_desc_nodes()
				# extra.gather_child_idx()

				# print parent_node.children[0].name
				# print extra.children[0].name
		else:
			extra = self.create_node('Extra', values=extra_idx)
			#extra cluster is the root
			extra.add_children(self.root)
			self.root = extra


def baseline_ad(scenario,print_ad=False, print_ccm = False):
	ad=get_ad_nvar(scenario,size_clusters=[3,2,2,3,2,4], extra_prop=2.0/12.0)
	ccm, clusters = get_ccm_nvar(scenario,size_clusters=[3,2,2,3,2,4], extra_prop=2.0/12.0)
	if print_ad == True:
		print ad
	if print_ccm == True:
		print ccm
	
	return ad,ccm

def make_truth():
	tree = Tree(name='N1', size=3)
	root = tree.root
	node2 = tree.create_node('N2',tree.root,2)
	node3 = tree.create_node('N3',tree.root,2)
	node4 = tree.create_node('N4',node2,3)
	node5 = tree.create_node('N5',node2,2)
	node6 = tree.create_node('N6',node3,4)
	return(tree)

def tree_from_df(tree_df, size = None):
	if size is None:
		size = np.repeat(2,tree_df.shape[0])
	tree = Tree(name = tree_df[0,0], size = size[0])
	nodes = dict()
	nodes[1] = tree.root		
	for i in range(1,tree_df.shape[0]):
		node = tree.create_node(tree_df[i,0],nodes[tree_df[i,1]], size = size[i])
		nodes[tree_df[i,0]] = node
	return(tree)


def tree_mistakes(tree):
	print "orig"
	print tree.tree_struct()
	# split_tree = split_bottom(tree)
	# merged_bottom_tree = merge_bottom(tree)
	# merged_top_tree = merge_top(tree)
	# extra_bottom_tree = add_extra_bottom(tree)
	# merged_extra_tree = merge_top_and_add_extra(tree)
	switched_tree = parent_is_grandparent(tree)
	# switched_tree = sibling_is_parent(tree)

def split_bottom(tree):
	tree_df = tree.tree_struct()
	split_bottom_tree = copy.deepcopy(tree)
	split_bottom_tree.split_node(max(tree_df[:,0]), max(tree_df[:,0])+1, same =False)
	split_node = split_bottom_tree.get_node(max(tree_df[:,0])+1)
	split_node.plot_name = max(tree_df[:,0])
	print "split"
	print split_bottom_tree.tree_struct(plot_labels = True)
	print "orig"
	print tree.tree_struct()
	return(split_bottom_tree)

def merge_bottom(tree):
	tree_df = tree.tree_struct()
	merged_bottom_tree = copy.deepcopy(tree)
	bottom_siblings = tree_df[tree_df[:,1] == tree_df[np.argmax(tree_df[:,0]),1],:]
	if bottom_siblings.shape[0] > 1:
		node_1 = bottom_siblings[0,0]
		node_2 = bottom_siblings[1,0]

	else:
		bottom_node = tree_df[np.argmax(tree_df[:,0]),:]
		print bottom_node
		node_1 = bottom_node[1]
		node_2 = bottom_node[0]
		
	merged_bottom_tree.merge_nodes(node_1, node_2)
	merged_node = merged_bottom_tree.get_node(node_1)
	merged_node.plot_name = str(node_1) +'/' + str(node_2)
	print merged_node.value
		
	print merged_bottom_tree.tree_struct(plot_labels = True)
	return(merged_bottom_tree)

def merge_top(tree):
	tree_df = tree.tree_struct()
	merged_top_tree = copy.deepcopy(tree)
	node_1 = 1
	node_2 = 2
		
	merged_top_tree.merge_nodes(node_1, node_2)
	merged_node = merged_top_tree.get_node(node_1)
	merged_node.plot_name = str(node_1) +'/' + str(node_2)
	print merged_node.value
		
	print merged_top_tree.tree_struct(plot_labels = True)
	return(merged_top_tree)

def add_extra_bottom(tree):
	tree_df = tree.tree_struct()
	parent_name = max(tree_df[:,0]) 
	extra_bottom_tree = copy.deepcopy(tree)
	extra_bottom_tree.extra_node(0.5,parent_name,all_nodes=True, num_nodes=None)
	parent_node = extra_bottom_tree.get_node('Extra')
	parent_node.plot_name = 'X'
	print extra_bottom_tree.tree_struct(plot_labels=True)
	return extra_bottom_tree

def merge_top_and_add_extra(tree):
	merged_tree = merge_top(tree)
	print tree.tree_struct()
	extra_and_merged = add_extra_bottom(merged_tree)

def parent_is_grandparent(tree):
	tree_df = tree.tree_struct()
	node_name = max(tree_df[:,0]) 
	pig_tree = copy.deepcopy(tree)
	node = pig_tree.get_node(node_name)
	grandparent = node.parent.parent
	if( len(grandparent.children) > 1):
		sibling_is_parent(tree)
	else:
		pig_tree.switch_parent(node_name, grandparent.name)
		print pig_tree.tree_struct(plot_labels = True)
		return pig_tree

def sibling_is_parent(tree):
	tree_df = tree.tree_struct()
	node_name = max(tree_df[:,0])
	sip_tree = copy.deepcopy(tree)
	node = sip_tree.get_node(node_name)
	print node.parent.name
	if (len(node.parent.children) > 1):

		sibling_name = tree_df[(tree_df[:,1] == node.parent.name) & (tree_df[:,0] != node_name),0]
	else:
		print "here"
		sibling_name = tree_df[(tree_df[:,1] == node.parent.parent.name) & (tree_df[:,0] != node.parent.name),0]
		node_name = node.parent.name
	print "sibling", sibling_name
	sip_tree.switch_parent(node_name, sibling_name)		
	print sip_tree.tree_struct(plot_labels = True)
	return sip_tree
	
def n_cluster_one_lineage(nssm, ordered =True):
	if ordered == True:
		tree = Tree(name ='N0', size =1)
		ssms = range(1,nssm)
	last_ssm  = tree.root
	for ssm in ssms:
		node = tree.create_node('N'+str(ssm), parent = last_ssm, size =1)
		last_ssm = node
	return tree

def one_cluster(nssm, ordered =True):
	tree = Tree(name ='N0', size = nssm)
	return tree

def n_cluster_two_lineages(nssm, ordered = True):
	if ordered == True:
		tree = Tree(name='N0', size =1)
		l1_root = tree.create_node('N1', parent = tree.root, size =1)
		l1_ssms = range(1,nssm/2)

	last_ssm  = l1_root
	for ssm in l1_ssms:
		node = tree.create_node('N'+str(ssm), parent = last_ssm, size =1)
		last_ssm = node

	l2_root = tree.create_node('N2', parent = tree.root, size =1)
	l2_ssms = range((nssm/2+1),(nssm-1))
	last_ssm  = l2_root
	for ssm in l2_ssms:
		node = tree.create_node('N'+str(ssm), parent = last_ssm, size =1)
		last_ssm = node
	return tree


def ncluster_correct_lineage(tree=None):
	if tree is None:
		tree = make_truth()
	cluster_last_node = {}
	cluster_tree = []
	for node_idx in range(len(tree.nodes)):
		node = tree.nodes[node_idx]
		last_node = None
		for idx in range(len(node.value)):
			parent = ''
			if node.parent is None and idx == 0:
				cluster_tree = Tree("N0",values = node.value[idx])
				last_node = cluster_tree.root
				continue
			else:
				if idx == 0:
					parent = cluster_last_node[node.parent.name]				
				else: 
					parent = last_node
			new_node = cluster_tree.create_node('N', parent = parent, values=node.value[idx])
			if idx == (len(node.value)-1):
				cluster_last_node[node.name] = new_node
			last_node = new_node
	# print cluster_tree.ccm()
	return cluster_tree
			
def split_cluster(tree = None, name_split=None, name_new=None, same=True):
	if tree is None:
		tree = make_truth()
	tree.split_node(name_split,name_new,same=same)
	n6 = tree.get_node(name_split)
	n7 = tree.get_node(name_new)
	return tree

def merge_clusters(node_name1=None, node_name2=None):
	tree = make_truth()
	tree.merge_nodes(node_name1, node_name2)
	n1 = tree.get_node(node_name1)
	return tree

def collapse_clusters(node_name):
	tree = make_truth()
	tree.collapse_node(node_name)
	return tree

def switch_parent_cluster(node_name=None, new_parent_name=None):
	tree = make_truth()
	tree.switch_parent(node_name, new_parent_name)
	#print tree.ad()
	return tree

def extra_cluster(extra_prop=1.0/6.0, parent_name=None, all_nodes = True, num_nodes = None):
	tree = make_truth()
	tree.extra_node(extra_prop,parent_name,all_nodes=all_nodes, num_nodes=num_nodes)
	return tree


########################New mistake scenarios


def branching_test():
	#truth2
	tree = Tree(name='N1', size=3)
	node2 = tree.create_node('N2',tree.root,2)
	node3 = tree.create_node('N3',tree.root,3)
	node4 = tree.create_node('N4',node3,3)
	return tree

def branching_test2():
	#truth2
	tree = Tree(name='N1', size=3)
	node2 = tree.create_node('N2',tree.root,2)
	node3 = tree.create_node('N3',tree.root,3)
	node4 = tree.create_node('N4',node2,3)
	node5 = tree.create_node('N5',node2,3)
	node6 = tree.create_node('N6',node3,3)
	node7 = tree.create_node('N7',node3,3)
	node8 = tree.create_node('N8',node6,3)
	return tree

def branching_test3():
	#truth2
	tree = Tree(name='N1', size=3)
	node2 = tree.create_node('N2',tree.root,2)
	node3 = tree.create_node('N3',tree.root,3)
	node4 = tree.create_node('N4',node2,3)
	node5 = tree.create_node('N5',node2,3)
	node6 = tree.create_node('N6',node3,3)
	node7 = tree.create_node('N7',node3,3)
	node8 = tree.create_node('N8',node4,3)
	return tree

def branching_test4():
	#truth2
	tree = Tree(name='N1', size=3)
	node3 = tree.create_node('N3',tree.root,3)
	node2 = tree.create_node('N2',tree.root,2)
	node4 = tree.create_node('N4',node2,3)
	node5 = tree.create_node('N5',node4,3)
	return tree

def branching_test5():
	#truth2
	tree = Tree(name='N1', size=3)
	node2 = tree.create_node('N2',tree.root,2)
	node3 = tree.create_node('N3',node2,3)
	node4 = tree.create_node('N4',tree.root,3)
	node5 = tree.create_node('N5',node4,3)
	node6 = tree.create_node('N6',node5,3)
	return tree

def branching_test6():
	#truth2
	tree = Tree(name='N1', size=3)
	node2 = tree.create_node('N2',tree.root,2)
	node3 = tree.create_node('N3',node2,3)
	node4 = tree.create_node('N4',node3,3)
	node5 = tree.create_node('N5',tree.root,3)
	return tree

def branching_test7():
	#truth2
	tree = Tree(name='N1', size=3)
	node2 = tree.create_node('N2',tree.root,2)
	node3 = tree.create_node('N3',node2,3)
	node4 = tree.create_node('N4',node3,3)
	node5 = tree.create_node('N5',tree.root,3)
	node6 = tree.create_node('N6',node2,3)
	return tree

def branching_test8():
	#truth2
	tree = Tree(name='N1', size=3)
	node2 = tree.create_node('N2',tree.root,2)
	node3 = tree.create_node('N3',node2,3)
	node4 = tree.create_node('N4',tree.root,3)
	node5 = tree.create_node('N5',node2,3)
	node6 = tree.create_node('N6',node3,3)
	return tree

def linear_to_branching(tree):
	#truth2 as a mistake for truth1
	tree.switch_parent('N4','N2')
	node = tree.get_node('N3')
	node.gather_child_idx()
	print tree.ad()

def matches(scenario, f, pass_l, scenarios,  **kwargs):
	orig_ad, orig_ccm = baseline_ad(scenario)
	new_tree = f(**kwargs)
	new_ad = new_tree.ad()

	new_ccm = new_tree.ccm()
	print scenario
	pass_v = 0
	if (np.array_equal(new_ad,orig_ad)):
	    print('Pass AD')
	    # return('Pass')
	else: 
		print('Fail AD')
		print "new"
		print new_ad
		print "orig"
		print orig_ad
		pass_v = 1
	if (np.array_equal(new_ccm,orig_ccm)):
	    print('Pass CCM')
	    # return('Pass')
	else: 
		print('Fail CCM')
		print "new"
		print new_ccm
		print "orig"
		print orig_ccm
		pass_v = 1
	pass_l.append((scenario, pass_v))



def test_all():
	pass_l = []
	scenarios = []
	matches('Truth',  make_truth, pass_l, scenarios)
	matches('SplitClusterBotSame', split_cluster, pass_l, scenarios, name_split='N6', name_new='N7',same=True)
	matches('SplitClusterBotDiff', split_cluster, pass_l, scenarios, name_split='N6', name_new='N7',same=False)
	matches('SplitClusterMidOneChild', split_cluster,  pass_l, scenarios,name_split='N3', name_new='N7',same=True)
	matches('SplitClusterMidMultiChild', split_cluster, pass_l, scenarios, name_split='N2', name_new='N7',same=True)
	matches('MergeClusterMid&BotOneChild', merge_clusters, pass_l, scenarios, node_name1='N3',node_name2='N6')
	matches('MergeClusterBot', merge_clusters, pass_l, scenarios, node_name1='N4',node_name2='N5')
	matches('MergeClusterTop&Mid', merge_clusters,  pass_l, scenarios,node_name1='N1',node_name2='N2')
	matches('MergeClusterMid&BotMultiChild', merge_clusters, pass_l, scenarios, node_name1='N2',node_name2='N5')
	matches('ParentIsSibling', switch_parent_cluster, pass_l, scenarios, node_name='N5',new_parent_name='N4')
	matches('ParentIsAunt', switch_parent_cluster, pass_l, scenarios, node_name='N5',new_parent_name='N3')
	matches('ParentIsCousin', switch_parent_cluster, pass_l, scenarios, node_name='N5',new_parent_name='N6')
	matches('ParentIsNieceWithChildren', switch_parent_cluster, pass_l, scenarios, node_name='N2',new_parent_name='N6')
	matches('ParentIsGrandparent', switch_parent_cluster, pass_l, scenarios, node_name='N5',new_parent_name='N1')
	matches('SmallExtraCurBot', extra_cluster, pass_l, scenarios, parent_name='N3')
	#original SmallExtraMid is wrong
	matches('SmallExtraMid', extra_cluster, pass_l, scenarios, parent_name='N1')
	matches('SmallExtraTop', extra_cluster, pass_l, scenarios, parent_name=None)
	matches('NClusterCorrectLineage', ncluster_correct_lineage, pass_l, scenarios)
	matches('NClusterOneLineage', n_cluster_one_lineage, pass_l, scenarios,nssm = 16)
	matches('NClusterTwoLineages', n_cluster_two_lineages, pass_l, scenarios,nssm=16)
	matches('OneCluster', one_cluster, pass_l, scenarios,nssm=16)
	fail = [sc[0] for sc in pass_l if sc[1] > 0]
	if len(fail)>0:
		for sc in fail:
			print "failed at", sc  #print scenarios#scenarios[np.where(pass_l>0)]
	else:
		print "all passed"

def test_scenario(t,f,**kwargs):
	tree = t()
	new_tree = f(tree, **kwargs)
	new_ad = new_tree.ad()
	new_ccm = new_tree.ccm()
	print new_ad
	print new_ccm


# print "linear truth"
# tree = linear_branching_truth()
# # tree = linear_truth()
# print tree.ad()
# print tree.ccm()
# print "mistake"
# tree = collapse_all_clusters_and_extra(tree)
# print tree.ad()
# print tree.ccm()

# print "branching truth"
# tree = branching_test8()
# tree.standard_node_naming()
# tree_df = tree.tree_struct(uniform=True, plot_labels = True, tier=True)

# n2 = tree.get_node('N2')
# n2.gather_desc_nodes()
# print n2.desc

# tree.switch_parent('N6','N3')
# tree.standard_node_naming()
# tree_df = tree.tree_struct(uniform=True, plot_labels = True, tier=True)

# n2 = tree.get_node('N2')
# n2.gather_desc_nodes()
# print n2.desc


# tree = branching_test8()
# tree.standard_node_naming()
# # tree.uniform_node_naming()
# tree_df = tree.tree_struct(uniform=True, plot_labels = True, tier=True)
# print tree_df
# # print tree.tiers
# # tree.gather_tiered_nodes()
# tier_three = tree.get_tier(3)



# map(lambda x: pprint.pprint(x.name) ,tier_three)
# pprint.pprint(tree.node_tiers)
# pprint.pprint(tier_three)

# tree.uniform_node_naming()
# tree_df = tree.tree_struct(plot_labels = True,uniform=True)
# print tree_df

# print tree.ad()
# print tree.ccm()
# print "mistake"
# tree = branching_to_linear(tree)
# print tree.ad()
# print tree.ccm()




# tree = make_truth()
# ncluster_correct_lineage()
# split(name_split='N2', name_new='N7',same=False)
# merge_clusters('N3','N6')
# baseline_ad('OneCluster', print_ad=False, print_ccm=True)
# tree = extra_cluster(parent_name=None)
# tree = switch_parent_cluster('N5','N6')
# print tree.ad()
# print tree.ccm()
# matches('NClusterCorrectLineage', ncluster_correct_lineage, [], [])
# tree = one_cluster(16)
# print tree.ad()
# print tree.ccm()
# baseline_ad('Truth', print_ccm=True)
# collapse_node('N2')
# test_all()



