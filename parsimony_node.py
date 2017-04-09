class Parsimony_Node:
	'''A node for an evolutionary tree, meant to be built from 
	newick format using character states from a nexus file'''

	def __init__(self,name,state = None):
		'''Node contains a character state, a name, 
		and pointers to left and right children'''
		self.state = state
		self.name = name
		self.right = None
		self.left = None

	def is_leaf(self):
		'''Returns true if the node is a leaf'''
		return self.left == None and self.right == None
