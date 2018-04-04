import numpy as np;
from make_ad_nvar import *
from make_ccm import *


class Node:
	def __init__(self,value, name, parent=None):
		self.value = value
		self.children = np.array([])
		self.desc = np.array([])
		self.parent = parent
		self.name = name
		self.desc_idx = np.array([],np.int32)
		if parent is not None:
			parent.add_children([self])
	def add_children(self, children):
		if self.children is not None:
			for child in children:
				self.children = np.append(self.children,child)
	def remove_child(self, child_name):
		for idx in range(len(self.children)):
			if self.children[idx].name == child_name:
				print "removing ", self.children[idx].name, idx
				self.children = np.delete(self.children, idx)
				#only want to remove one child at a time or the loop breaks
				return None
	def gather_desc_nodes(self, top_node= None):
		if top_node is None:
			top_node = self
			for child in self.children:
				top_node.desc = np.append(top_node.desc, child.name)
				if len(child.children) > 0:
					#get the indices for the descendants of each child
					child.gather_desc_nodes(top_node)	
	def gather_child_idx(self,top_node=None):
		#top node is the node for which we will gather all of the descendant indices
		if top_node is None:
			top_node = self
		for child in self.children:
			#get the indices for each child
			top_node.desc_idx = np.concatenate((top_node.desc_idx, child.value))
			if len(child.children) > 0:
				#get the indices for the descendants of each child
				child.gather_child_idx(top_node)
	def collapse_children(self, full=True):
		#ALL of the node's descendants (not just children) are merged into the cluster
		if full == True:
			top_node = self
			self.gather_child_idx(top_node)
			#all of the node's descendants now become part of the cluster
			self.value = np.append(self.value,self.desc_idx)
			self.desc_idx = np.array([],np.int32)
			self.children = np.array([])

class Tree:
	def __init__(self,size,name):
		self.nodes = []
		root_idx = self.node_indices(size)
		self.root = Node(root_idx,name)
		self.nodes = [self.root]
	def node_indices(self,size):
		idx = np.arange(size)
		nssm = self.comp_nssms()
		idx = idx + nssm
		return idx
	def create_node(self,name,parent=None,size=None, values=None):
		if size is not None and values is None:
			idx = self.node_indices(size)
		if values is not None:
			idx = values
		node = Node(idx, name, parent)
		self.nodes.append(node)
		return node
	def comp_nssms(self):
		nssms = reduce(lambda x,y: x + len(y.value),self.nodes ,0)
		return(nssms)
	def ccm(self):
		nssms = self.comp_nssms()
		out = np.zeros((nssms,len(self.nodes)))
		print nssms
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
			node.gather_child_idx()
			desc_idx = node.desc_idx
			out[np.ix_(node.value, desc_idx)] = 1
		return out
	def get_node(self,name, return_idx=False):
		print len(self.nodes)
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

	def merge_nodes(self, node_name1, node_name2, collapse = False):
		node1 = self.get_node(node_name1)
		node2, idx2 = self.get_node(node_name2, return_idx=True)
		if collapse == False:
			node1.value = np.concatenate((node1.value,node2.value))
			node1.remove_child(node_name2)
			if len(node2.children) > 0 :
				node1.add_children(node2.children)
		elif collapse == True:
			print "before", node2.value
			node2.collapse_children()
			print "after", node2.value
			node1.value = np.concatenate((node1.value,node2.value))
			print "merged", node1.value
		self.nodes = np.delete(self.nodes, idx2)
		# print self.nodes
		# print node1.values 
	def collapse_node(self, node_name):
		node = self.get_node(node_name)
		node.gather_desc_nodes()
		node.collapse_children()
		for name in node.desc:
			d_node, d_idx = self.get_node(name,return_idx = True)
			self.nodes = np.delete(self.nodes, d_idx)

def baseline_ad(scenario,print_ad=False, print_ccm = False):
	ad=get_ad_nvar(scenario,size_clusters=[3,2,2,3,2,4])
	ccm, clusters = get_ccm_nvar(scenario,size_clusters=[3,2,2,3,2,4])
	if print_ad == True:
		print ad
	if print_ccm == True:
		print ccm
	
	return ad,ccm

def make_truth():
	tree = Tree(3,'N1')
	root = tree.root
	node2 = tree.create_node('N2',tree.root,2)
	node3 = tree.create_node('N3',tree.root,2)
	node4 = tree.create_node('N4',node2,3)
	node5 = tree.create_node('N5',node2,2)
	node6 = tree.create_node('N6',node3,4)
	return(tree)

def split_cluster(name_split=None, name_new=None, same=True):
	tree = make_truth()
	print "same", same
	tree.split_node(name_split,name_new,same=same)
	n6 = tree.get_node(name_split)
	print n6.value
	n7 = tree.get_node(name_new)
	print n7.value
	print len(tree.nodes)
	ad = tree.ad()

	# print tree.ccm()
	print ad
	return tree
def merge_clusters(node_name1=None, node_name2=None):
	tree = make_truth()
	print(len(tree.nodes))
	tree.merge_nodes(node_name1, node_name2)
	n1 = tree.get_node(node_name1)
	print n1.value
	print n1.name
	print n1
	#print tree.nodes
	print tree.ccm()
	# print tree.ad()
	return tree

def collapse_node(node_name):
	tree = make_truth()
	tree.collapse_node(node_name)
	print tree.ccm()
	print tree.ad()
	return tree

def matches(scenario, f, **kwargs):
	orig_ad, orig_ccm = baseline_ad(scenario)
	new_tree = f(**kwargs)
	new_ad = new_tree.ad()

	new_ccm = new_tree.ccm()
	if (np.array_equal(new_ad,orig_ad)):
	    print('Pass AD')
	    # return('Pass')
	else: 
		print('Fail AD')
		print new_ad
		print orig_ad
	if (np.array_equal(new_ccm,orig_ccm)):
	    print('Pass CCM')
	    # return('Pass')
	else: 
		print('Fail CCM')
		print new_ccm
		print orig_ccm

# split(name_split='N2', name_new='N7',same=False)
# merge_clusters('N3','N6')
# baseline_ad('MergeClusterMid&BotMultiChild', print_ccm=True)
# baseline_ad('Truth', print_ccm=True)
# collapse_node('N2')
# matches('Truth', make_truth)
#matches('SplitClusterBotDiff', split_cluster_bot, same=False)
# matches('SplitClusterBotSame', split_cluster_bot, name_split='N6', name_new='N7',same=True)
# matches('SplitClusterBotDiff', split_cluster_bot, name_split='N6', name_new='N7',same=False)
#matches('SplitClusterMidOneChild', split_cluster_bot, name_split='N3', name_new='N7',same=True)
#matches('SplitClusterMidMultiChild', split_cluster, name_split='N2', name_new='N7',same=True)
# matches('MergeClusterMid&BotOneChild', merge_clusters, node_name1='N3',node_name2='N6')
# matches('MergeClusterBot', merge_clusters, node_name1='N4',node_name2='N5')
# matches('MergeClusterTop&Mid', merge_clusters, node_name1='N1',node_name2='N2')
matches('MergeClusterMid&BotMultiChild', merge_clusters, node_name1='N2',node_name2='N5')