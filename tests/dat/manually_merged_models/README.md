This folder contains examples of merged models.

Each example has its own folder, structured as follows:
- M1_***.xml --> an SBML file of a model to be merged
- M2_***.xml --> an SBML file of the other model to be merged. It has a different namespace w.r.t. M1
- Merged_***.xml --> the resulting merged model
- Loopless_***.xml --> the resulting merged model after removing all ATP-generating cycles if there were some in the merged model
- met_matches.csv --> a csv file listing the matching metabolites in their 2 different namespaces. The column containing the metabolite ids under the namespace of model M1 is called M1_[name_of_namespace], same applies for M2. A third column, called ProtonsToAddToM2, indicates how many protons needs to be added to a molecule in M2's namespace, to correspond to the same metabolite under M1's namespace. An optional fourth column would be ProtonsToAddToM1 - same as the one before just for M1. 
- README.txt --> states where the models come from and who has merged them
