from Bio import AlignIO
from parsimony_tree import Parsimony_Tree
import numpy as np

class Parsimony_Scorer:
	'''Takes in a file of newick-format tree-strings
		and a nexus file of characters. Scores the newick
		trees against the entire dataset and finds the 
		most parsimonious tree.'''

	def __init__(self, newick_file, nexus_file):
		'''Opens the newick and nexus files and prepares data for analysis'''
		newick_file = open(newick_file,'rU')
		#assumes the newick file is formatted with each tree on its own line
		self.newick_trees = newick_file.readlines()

		self.nexus_file_name = nexus_file
		nexus_file = AlignIO.read(open("nexus_files/"+nexus_file),
			"nexus")
		self.nexus_matrix, self.nexus_taxa_order = self.set_chars_and_taxa_order(nexus_file)

		self.scored_trees = {}

	def set_chars_and_taxa_order(self,nexus_file):
		'''Sets the character matrix and taxa order from the nexus file'''
		nex_matrix = []
		nex_taxa = []
		for record in nexus_file:
			sequence = record.seq
			nex_matrix.append(list(str(sequence)))
			nex_taxa.append(record.name)
		return np.array(nex_matrix),nex_taxa

	def run(self):
		'''Runs the scorer over the data'''
		print("Finding best tree for " + self.nexus_file_name + "...")
		self.score_all_trees()
		self.find_best_tree()

	def score_all_trees(self):
		'''Scores all the trees against the dataset and collects the data
		in the dictionary self.scored_trees. The dictionary is keyed by score'''
		for newick in self.newick_trees:
			score = self.score_tree(newick)
			self.scored_trees[score] = newick

	def score_tree(self, newick):
		'''Scores a single tree against the entire dataset'''
		tree = Parsimony_Tree(newick)
		newick_taxa = tree.get_taxa_list()

		num_cols = len(self.nexus_matrix[0])
		total_parsimony_score = 0
		for i in range(num_cols):
			char_dict = self.get_character_dict(newick_taxa, i)
			tree.add_leaf_states(char_dict)
			node_list = tree.get_post_order_nodes()
			total_parsimony_score += self.fitch_bottom_up(node_list)
		return total_parsimony_score

	def get_character_dict(self, newick_taxa, i):
		'''Returns a dictionary of taxa-name and taxa-character
		to be sent to the tree's add_leaf_states method.
		Adapted from code by Carolyn Sy'''
		chars = self.nexus_matrix[:,i]
		char_dict = dict(zip(self.nexus_taxa_order, chars))
		return char_dict

	def fitch_bottom_up(self, node_list):
		'''Scores the tree using the bottom-up half of Fitch's Algorithm'''
		parsimony_score = 0
		for node in node_list:
			if not node.is_leaf():
				parsimony_score += self.set_internal_state(node)
		return parsimony_score

	def set_internal_state(self,node):
		'''Sets the state options for internal nodes'''
		l_state = node.left.state
		r_state = node.right.state
		intersection = self.get_intersection(l_state, r_state)
		if intersection:
			node.state = intersection
			return 0
		union = self.get_union(l_state, r_state)
		node.state = union
		return 1
	
	def get_intersection(self, l_state,r_state):
		'''takes in two node-states (either strings or lists)
		returns a list of the intersection of the states, or the
		state itself if there is only one item'''
		intersection = list(set(l_state).intersection(set(r_state)))
		if len(intersection) == 1:
			return intersection[0]
		return intersection

	def get_union(self, l_state,r_state):
		'''takes in two node-states (either strings or lists)
		returns a list of the union of the states, or the
		state itself if there is only one item'''
		union = list(set(l_state).union(set(r_state)))
		if len(union) == 1:
			return union[0]
		return union

	def find_best_tree(self):
		'''Finds and prints the best score and best scored tree.
		NOTE: if there are multiple trees with the same best score,
		you may not get the same tree returned each time.'''
		best_score = min(self.scored_trees.keys())
		print(best_score)
		print(self.scored_trees[best_score])

def main():
	scorer = Parsimony_Scorer("105treesOutgroupNamed.txt","morph_data.nex")
	scorer.run()
	scorer = Parsimony_Scorer("105treesOutgroupNamed.txt","RAG1_trimmed.nex")
	scorer.run()
	scorer = Parsimony_Scorer("105treesOutgroupNamed.txt","CYTB_trimmed.nex")
	scorer.run()


if __name__ == "__main__":
	main()
