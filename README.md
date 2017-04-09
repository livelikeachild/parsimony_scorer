# parsimony_scorer

`Parsimony_Scorer` takes in a file of newick-format strings
and a nexus file of data. It then builds and scores the given evolutionary trees, finding the most parsimonious tree according to the given data set.

The test data and function runs all possible rooted trees with 6 leaves on 3 datasets from separate nexus files.

### NOTES:

- `Parsimony_Tree` expects your newick strings to have spaces after commas. You have to manually update what's inside the `.split()` method on **line 32** if you want to remove them
- `Parsimony_Scorer` expects your newick file to have each tree on its own line
- Put any nexus files you use in the `nexus_files` folder, or manually remove `"nexus_files/"+` from **line 18** of `parsimony_scorer.py`
-  The names you use in your newick file must _match_ the names you use in your nexus files. 
- If more than one tree shares the best score, the code will not always return the same one. This is because the score-tree pairs are stored in a dictionary (`scored_trees`), which is ordered differently each time. This is why the score prints out _with_ the best tree.
- `rooted_named_trees.txt` is for testing the `Parsimony_Scorer` used in congruence with the nexus files in the `nexus_files` folder