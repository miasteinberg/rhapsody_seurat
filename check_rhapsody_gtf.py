import sys
import gffutils
from collections import defaultdict

# This script is just a collection of known conditions that will prevent a gtf
# file from running on the Rhapsody pipeline.
#
# It does not fix the problem, just identifies them.

# It accepts a gtf file and will print the results to std out
gtf = sys.argv[1]
output_list = []



#### Various problems and conditions to check for ####

# PROBLEM - gene_ids cannot have more than one gene_name associated with it
# SOLUTION - create a lookup dictionary to track gene_names
gene_lookup = defaultdict(list)

# PROBLEM - all features must have a gene_name
def check_name_exist(feature):
    if 'gene_name' not in feature.attributes:
        return output_list.append("One or more features don't have a gene_name")

# PROBLEM - all features must have strand information
def check_strand_exist(feature):
    if feature.strand == '.':
        return output_list.append("One or more features don't have strand information")




#### Iterate through the gtf using gffutils ####
# attributes for each feature are stored as lists
for f in gffutils.DataIterator(gtf):
    # features must have a gene_name
    check_name_exist(f)

    # features must have strand information
    check_strand_exist(f)

    # add to the gene_lookup dictionary
    gene_id = f.attributes['gene_id'][0]
    if 'gene_name' in f.attributes:
        name = f.attributes['gene_name'][0]
        gene_lookup[gene_id].append(name)


# check the lookup dictionary for gene_ids with more than one name
for value in gene_lookup.values():
    if len(list(set(value))) > 1:
        output_list.append('One or more gene_ids have more than one gene_name')



#### Print any problems that were found ####
if len(output_list) > 0:
    print('\n'.join([x for x in list(set(output_list))]))
else:
    print('No problems found!!')
