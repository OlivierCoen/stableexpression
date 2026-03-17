import pandas as pd
from genorm import m_measure, genorm
import sys

file = sys.argv[1]
counts = pd.read_csv(file, header=0, decimal=",")
counts = counts.drop(columns=["groupe"]).set_index("echantillon")
print(counts)

# Compute `M` value for this set of control genes.
m_measure(counts)

# Select top 2 control genes with lowest `M`.
gene_names, m_values = genorm(counts, n_stop=11)
print(gene_names, m_values)
