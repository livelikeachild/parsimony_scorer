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
		self._newick_trees = newick_file.readlines()

		self._nexus_file_name = nexus_file
		nexus_file = AlignIO.read(open("nexus_files/"+nexus_file),
			"nexus")
		self._nexus_matrix, self._nexus_taxa_order = self._set_chars_and_taxa_order(nexus_file)

		self._scored_trees = {}

	def _set_chars_and_taxa_order(self,nexus_file):
		'''Sets the character matrix and taxa order from the nexus file'''
		nex_matrix = []
		nex_taxa = []
		for record in nexus_file:
			sequence = record.seq
			nex_matrix.append(list(str(sequence)))
			nex_taxa.append(record.name)
		return np.array(nex_matrix),nex_taxa

	def update_taxa_names(self,taxa_names):
		'''If your newick-file taxa names don't match those in the nexus
		file, give the newick-file names here'''
		self._nexus_taxa_order = taxa_names

	def get_best_tree(self):
		'''Runs the scorer over the data'''
		self._score_all_trees()
		best_tree = self._find_best_tree()
		return best_tree

	def _score_all_trees(self):
		'''Scores all the trees against the dataset and collects the data
		in the dictionary self.scored_trees. The dictionary is keyed by score'''
		for newick in self._newick_trees:
			score = self._score_tree(newick)
			self._scored_trees[score] = newick

	def _score_tree(self, newick):
		'''Scores a single tree against the entire dataset'''
		tree = Parsimony_Tree(newick)

		num_cols = len(self._nexus_matrix[0])
		total_parsimony_score = 0
		for column in range(num_cols):
			char_dict = self._get_character_dict(column)
			tree.add_leaf_states(char_dict)
			node_list = tree.get_post_order_nodes()
			total_parsimony_score += self._fitch_bottom_up(node_list)
		return total_parsimony_score

	def _get_character_dict(self, column):
		'''Returns a dictionary of taxa-name and taxa-character
		to be sent to the tree's add_leaf_states method.
		Adapted from code by Carolyn Sy'''
		chars = self._nexus_matrix[:,column]
		char_dict = dict(zip(self._nexus_taxa_order, chars))
		return char_dict

	def _fitch_bottom_up(self, node_list):
		'''Scores the tree using the bottom-up half of Fitch's Algorithm'''
		parsimony_score = 0
		for node in node_list:
			if not node.is_leaf():
				parsimony_score += self._find_internal_state(node)
		return parsimony_score

	def _find_internal_state(self,node):
		'''Sets the state options for internal nodes'''
		l_state = node.get_left().get_state()
		r_state = node.get_right().get_state()
		intersection = self._get_intersection(l_state, r_state)
		if intersection:
			node.set_state(intersection)
			return 0
		union = self._get_union(l_state, r_state)
		node.set_state(union) 
		return 1
	
	def _get_intersection(self, l_state,r_state):
		'''takes in two node-states (either strings or lists)
		returns a list of the intersection of the states, or the
		state itself if there is only one item'''
		intersection = list(set(l_state).intersection(set(r_state)))
		if len(intersection) == 1:
			return intersection[0]
		return intersection

	def _get_union(self, l_state,r_state):
		'''takes in two node-states (either strings or lists)
		returns a list of the union of the states, or the
		state itself if there is only one item'''
		union = list(set(l_state).union(set(r_state)))
		if len(union) == 1:
			return union[0]
		return union

	def _find_best_tree(self):
		'''Finds and prints the best score and best scored tree.
		NOTE: if there are multiple trees with the same best score,
		you may not get the same tree returned each time.'''
		best_score = min(self._scored_trees.keys())
		print("Best Parsimony score:",best_score)
		return self._scored_trees[best_score]

def main():
	nexus_files = ["morph_data.nex","RAG1_trimmed.nex","CYTB_trimmed.nex"]

	for nexus in nexus_files:
		scorer = Parsimony_Scorer("rooted_trees.txt",nexus)
		scorer.update_taxa_names(['0','1','2','3','4','5'])
		print("Finding best tree for " + nexus + "...")
		best_tree = scorer.get_best_tree()
		print("Best Tree: ",best_tree)

	for nexus in nexus_files:
		scorer = Parsimony_Scorer("rooted_named_trees.txt",nexus)
		print("Finding best tree for " + nexus + "...")
		best_tree = scorer.get_best_tree()
		print("Best Tree: ",best_tree)

if __name__ == "__main__":
	main()
