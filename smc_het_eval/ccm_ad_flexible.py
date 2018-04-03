import numpy as np;

class Node:
	def __init__(self,value , name, parent=None):
		self.value = value
		self.children = []
		self.parent = parent
		self.name = name
		if parent is not None:
			parent.add_children([self])
	#	else:
	#		self.value = [0:size]
	def add_children(self, children):
		if self.children is not None:
			for child in children:
				self.children.append(child)

class Tree:
	def __init__(self,value,size,name):
		self.nodes = []
		root_idx = self.node_indices(size)
		self.root = Node(root_idx,'root')
		self.nodes = [self.root]
	def node_indices(self,size):
		idx = np.arange(size)
		nssm = self.comp_nssms()
		idx = idx + nssm
		print "idx",idx
		return idx
	def create_node(self,value,size,name,parent):
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
		print out
		ccm = np.dot(out,out.T)
		return ccm




tree = Tree([0,1,2],3,'root')
root = tree.root
child_1 = tree.create_node([3,4],2,'child1',tree.root)
#root = Tree()
#print(root.value)
print "child1 value ", child_1.value
print(child_1.parent)
print(root.children)
print len(child_1.value)
print root is child_1.parent
print tree.ccm()