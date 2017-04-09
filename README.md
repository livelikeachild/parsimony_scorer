# parsimony_scorer
Builds and scores evolutionary trees based on newick-format strings and nexus files of data

Parsimony_Scorer takes in a file of newick-format strings (expects one tree on each line)
and a nexus file of data. 

The test data and function runs all possible rooted trees with 6 leaves on 3 datasets from separate nexus files.

###NOTES:

- `parsimony_tree` expects your newick strings to have spaces after commas. You have to manually update what's inside the `.split()` method on line 32 if you want to remove them
-  The names you use in your newick file must _match_ the names you use in your nexus files. 
- If more than one tree shares the best score, the code will not always return the same one. This is because the score-tree pairs are stored in a dictionary (`scored_trees`), which is ordered differently each time. This is why the score prints out _with_ the best tree. 