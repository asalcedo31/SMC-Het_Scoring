import numpy as np;
from make_ad_nvar import *
from make_ccm import *


class Node:
	def __init__(self,value, name, parent=None):
		self.value = value
		self.children = []
		self.parent = parent
		self.name = name
		self.desc_idx = np.array([],np.int32)
		if parent is not None:
			parent.add_children([self])
	def add_children(self, children):
		if self.children is not None:
			for child in children:
				self.children.append(child)
	def gather_child_idx(self,top_node=None):
		if top_node is None:
			top_node = self
		for child in self.children:
			top_node.desc_idx = np.concatenate((top_node.desc_idx, child.value))
			if len(child.children) > 0:
				child.gather_child_idx(top_node)

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
	def create_node(self,size,name,parent):
		idx = self.node_indices(size)
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
		ccm = np.dot(out,out.T)
		return ccm
	def gather_children(self):
		desc_idx = []

	def ad(self):
		nssms = self.comp_nssms()
		out = np.zeros((nssms,nssms))
		for node in self.nodes:
			node.gather_child_idx()
			desc_idx = node.desc_idx
			out[np.ix_(node.value, desc_idx)] = 1
		return out



def baseline_ad(scenario):
	ad=get_ad_nvar(scenario,size_clusters=[3,2,2,3,2,2])
	ccm, clusters = get_ccm_nvar(scenario,size_clusters=[3,2,2,3,2,2])
	print ccm
	return ad,ccm

def make_truth():
	tree = Tree(3,'N1')
	root = tree.root
	node2 = tree.create_node(2,'N2',tree.root)
	node3 = tree.create_node(2,'N3',tree.root)
	node4 = tree.create_node(3,'N4',node2)
	node5 = tree.create_node(2,'N5',node2)
	node6 = tree.create_node(2,'N6',node3)
	print tree.ccm()
	return(tree.ad(), tree.ccm())

def matches(scenario, f):
	orig_ad, orig_ccm = baseline_ad(scenario)
	new_ad, new_ccm = f()
	if (np.array_equal(new_ad,orig_ad)):
	    print('Pass AD')
	    # return('Pass')
	else: 
	    print('Fail AD')
	if (np.array_equal(new_ccm,orig_ccm)):
	    print('Pass CCM')
	    # return('Pass')
	else: 
	    print('Fail CCM')
matches('Truth', make_truth)