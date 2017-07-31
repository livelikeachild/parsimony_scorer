from parsimony_node import Parsimony_Node

class Parsimony_Tree:

	def __init__(self, newick, internal = 7):
		'''Initializes an evolutionary tree and populates 
		it based on the given newick string with optional 
		"internal" value to name the internal nodes. If no 
		number is given, it assumes that the tree will
		have 6 leaves (with internal starting at 7)'''
		self._root = None
		self._newick = newick
		self._internal = internal
		self._leaf_num = 0
		self._leaves= []
		self._taxa_list = self._set_taxa_list()

		self._populate_parsimony_tree()

	def _set_taxa_list(self):
		'''Sets the taxa list to hold the names of the taxa in the 
		order they appear in the newick string'''
		removals = ['(',')',',',';']
		taxa_string = self._newick
		for char in removals:
			taxa_string = taxa_string.replace(char,"")

		return taxa_string.split()

	def _populate_parsimony_tree(self):
		'''Builds the tree from the newick string. 
		Assumes your newick string contains a space after each 
		comma and splits it accordingly. Change what's in split's 
		()'s below to match your newick string.'''
		newick = self._newick.split(', ')
		self._root, newick = self._add_node(newick)

	def _add_node(self, newick):
		'''Recursive function that builds an entire tree from a 
		newick string. If chars data was given, it will also set 
		the leaf character states'''
		
		if newick[0][0] == "(":
			#must create an internal node
			node = Parsimony_Node(self._internal)
			self._internal += 1
			
			#get rid of the '(' and create the internal node's left child
			newick = [newick[0][1:]]+newick[1:]
			node.left, newick = self._add_node(newick)
			#get rid of the left child and create the right
			node.right, newick = self._add_node(newick[1:])
		else:
			node = self._create_leaf_node()

		return node, newick

	def _create_leaf_node(self):
		'''Builds a leaf node using the taxa-list and characters
		if they exist'''
		node = Parsimony_Node(self._taxa_list[self._leaf_num])
		self._leaf_num += 1	
		self._leaves.append(node)
		return node

	def add_leaf_states(self,chars):
		'''The same tree will need to be scored for multiple characters.
		This resets the leaves' character states from a new dict of data.'''
		for node in self._leaves:
			node.state = chars[node.name]

	def get_post_order_nodes(self):
		'''Returns list of nodes in the tree in post-order traversal.'''
		self._post_order_nodes = []
		self._post_order_traversal_2(self._root)
		return self._post_order_nodes

	def _post_order_traversal_2(self,current_node):
		'''continues until the current node is null'''
		if current_node:
			self._post_order_traversal_2(current_node.left)
			self._post_order_traversal_2(current_node.right)
			self._post_order_nodes.append(current_node)
