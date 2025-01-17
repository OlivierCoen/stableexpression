import pandas as pd
import numpy as np


# Expression data for three control genes.
counts = pd.read_csv(
    "/home/olivier/repositories/nf-core-stableexpression/tests/input/pairwise_gene_variation/counts.csv",
    header=0,
)
counts.to_parquet(
    "/home/olivier/repositories/nf-core-stableexpression/tests/input/pairwise_gene_variation/counts.parquet",
    index=False,
)
counts.set_index("ensembl_gene_id", inplace=True)
counts = counts.T
print(counts)


def _m_numpy(gene_expression: np.ndarray) -> np.ndarray:
    """Internal control gene-stability measure `M`.

    Computes Eq. (4) in Ref. [1].

    [1]: Vandesompele, Jo, et al. "Accurate normalization of real-time quantitative
    RT-PCR data by geometric averaging of multiple internal control genes." Genome
    biology 3.7 (2002): 1-12.
    """

    if not (gene_expression > 0).all():
        raise ValueError(
            "Expression domain error: not all expression data are strictly positive!"
        )

    a = gene_expression
    # Eq. (2): A_{jk}^{(i)} = log_2 (a_{ij} / a_{ik})
    A = np.log2(np.einsum("ij,ik->ijk", a, 1 / a))
    print()
    print(A)
    print()
    # Eq. (3)
    V = np.std(A, axis=0)
    # Eq. (4) N.B., Since V_{j=k} is zero, we can simply ignore it since it does not
    # contribute to calculation.
    n = V.shape[1]
    return np.sum(V, axis=1) / (n - 1)


def m_measure(gene_expression):
    """Internal control gene-stability measure `M` as described in Ref. [1].

    [1]: Vandesompele, Jo, et al. "Accurate normalization of real-time quantitative
    RT-PCR data by geometric averaging of multiple internal control genes." Genome
    biology 3.7 (2002): 1-12.

    Args:
        gene_expression: Gene expression counts of `m` samples (rows) and `n` internal
            control genes (columns). Expression must be strictly positive.

    Raises:
        ValueError: Expression not strictly positive.
    """
    if isinstance(gene_expression, pd.DataFrame):
        m_values = _m_numpy(gene_expression.to_numpy())
        return pd.Series(m_values, index=gene_expression.columns)
    return _m_numpy(gene_expression)


print(m_measure(counts))
