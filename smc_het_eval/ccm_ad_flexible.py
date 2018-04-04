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
		if isinstance(child_names, basestring):
			child_names = [child_names]
		#list of indices to remove
		to_delete = np.array([],np.int32)
		for idx in range(len(self.children)):
			if self.children[idx].name in child_names:
				print "removing ", self.children[idx].name, idx
				to_delete = np.append(to_delete, idx)
		self.children = np.delete(self.children, to_delete)
		return None
	def gather_desc_nodes(self, top_node= None):
		if top_node is None:
			top_node = self
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
		node.gather_desc_nodes()
		node.collapse_children()
		for name in node.desc:
			d_node, d_idx = self.get_node(name,return_idx = True)
			self.nodes = np.delete(self.nodes, d_idx)

	def switch_parent(self, node_name,new_parent_name):
		node = self.get_node(node_name)
		print "old parent",node.parent.name
		print "new parent", new_parent_name
		new_parent_node = self.get_node(new_parent_name)
		old_parent_node = self.get_node(node.parent.name)
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

	def extra_node(self,extra_prop, parent_name, all_nodes=True, num_nodes = 2):
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
			# print node.name, node.value, num_taken
			# print taken_idx, taken
			extra_idx = np.append(extra_idx, taken)
			node.value = np.delete(node.value,taken_idx)
		# 	print "after", node.value
		# # print "extra", extra_idx
		# print self.nodes[1].children[0].name
		# print self.nodes[1].children[0].value
		# self.nodes[1].gather_child_idx()
		# print self.nodes[1].desc_idx
		if parent_name is not None:
			parent_node = self.get_node(parent_name)
			self.create_node('Extra', parent=parent_node,values=extra_idx)
		else:
			extra = self.create_node('Extra', values=extra_idx)
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
	tree.split_node(name_split,name_new,same=same)
	n6 = tree.get_node(name_split)
	n7 = tree.get_node(name_new)
	ad = tree.ad()
	return tree

def merge_clusters(node_name1=None, node_name2=None):
	tree = make_truth()
	tree.merge_nodes(node_name1, node_name2)
	n1 = tree.get_node(node_name1)
	return tree

def collapse_node(node_name):
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
	pass_l.append(pass_v)
	scenarios.append(scenario)



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



	print pass_l, scenarios
	if any(pass_l)>0:
		print "failed at", scenarios[pass_l>0]
	else:
		print "all passed"
# split(name_split='N2', name_new='N7',same=False)
# merge_clusters('N3','N6')
# baseline_ad('SmallExtraTop', print_ad=True, print_ccm=False)
# tree = extra_cluster(parent_name=None)
# tree = switch_parent_cluster('N5','N6')
# print tree.ad()
# print tree.ccm()
matches('SmallExtraTop', extra_cluster, [], [], parent_name=None)

# baseline_ad('Truth', print_ccm=True)
# collapse_node('N2')
# test_all()


