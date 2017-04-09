from parsimony_node import Parsimony_Node

class Parsimony_Tree:

	def __init__(self, newick, internal = 7):
		'''Initializes an evolutionary tree and populates 
		it based on the given newick string with optional 
		"internal" value to name the internal nodes. If no 
		number is given, it assumes that the tree will
		have 6 leaves (with internal starting at 7)'''
		self.__root = None
		self.__newick = newick
		self.__internal = internal
		self.__leaf_num = 0
		self.__leaves= []
		self.__taxa_list = self.__set_taxa_list()

		self.__populate_parsimony_tree()

	def __set_taxa_list(self):
		'''Sets the taxa list to hold the names of the taxa in the 
		order they appear in the newick string'''
		removals = ['(',')',',',';']
		taxa_string = self.__newick
		for char in removals:
			taxa_string = taxa_string.replace(char,"")

		return taxa_string.split()

	def __populate_parsimony_tree(self):
		'''Builds the tree from the newick string. 
		Assumes your newick string contains a space after each 
		comma and splits it accordingly. Change what's in split's 
		()'s below to match your newick string.'''
		newick = self.__newick.split(', ')
		self.__root, newick = self.__add_node(newick)

	def __add_node(self, newick):
		'''Recursive function that builds an entire tree from a 
		newick string. If chars data was given, it will also set 
		the leaf character states'''
		
		if newick[0][0] == "(":
			#must create an internal node
			node = Parsimony_Node(self.__internal)
			self.__internal += 1
			
			#get rid of the '(' and create the internal node's left child
			newick = [newick[0][1:]]+newick[1:]
			node.left, newick = self.__add_node(newick)
			#get rid of the left child and create the right
			node.right, newick = self.__add_node(newick[1:])
		else:
			node = self.__create_leaf_node()

		return node, newick

	def __create_leaf_node(self):
		'''Builds a leaf node using the taxa-list and characters
		if they exist'''
		node = Parsimony_Node(self.__taxa_list[self.__leaf_num])
		self.__leaf_num += 1	
		self.__leaves.append(node)
		return node

	def add_leaf_states(self,chars):
		'''The same tree will need to be scored for multiple characters.
		This resets the leaves' character states from a new dict of data.'''
		for node in self.__leaves:
			node.state = chars[node.name]

	def get_post_order_nodes(self):
		'''Returns list of nodes in the tree in post-order traversal.'''
		self.__post_order_nodes = []
		self.__post_order_traversal_2(self.__root)
		return self.__post_order_nodes

	def __post_order_traversal_2(self,current_node):
		'''continues until the current node is null'''
		if current_node:
			self.__post_order_traversal_2(current_node.left)
			self.__post_order_traversal_2(current_node.right)
			self.__post_order_nodes.append(current_node)
