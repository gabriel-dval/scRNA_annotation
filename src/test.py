'''
Utility file to make some helpful functions across all annotation
methods.
'''

import anndata
import pandas as pd
from scipy.io import mmread

def convert_mtx_to_h5ad(mtx, genes, barcodes, output_file):
    '''
    Convert a sparse matrix to a h5ad file.

    Parameters
    ----------
    mtx : scipy.sparse.csr_matrix
        The sparse matrix to convert.
    genes : list of str
        The gene names.
    barcodes : list of str
        The barcode names.
    output_file : str
        The output file path.

    Returns
    -------
    None
    '''
    # Load the matrix (MTX format)
    matrix = mmread(mtx).T.tocsr()

    # Load barcodes and genes
    barcodes = pd.read_csv(barcodes, header=None, sep="\t")[0].values
    genes = pd.read_csv(genes, header=None, sep="\t")[0].values

    adata = anndata.AnnData(matrix, obs=barcodes, var=genes)
    adata.obs.columns = adata.obs.columns.astype(str)
    adata.var.columns = adata.var.columns.astype(str)
    adata.write(output_file)


if __name__ == '__main__':
    mtx = '../data/raw_sbm/raw_sbm_data_2024.11.15/matrix.mtx'
    genes = '../data/raw_sbm/raw_sbm_data_2024.11.15/genes.tsv'
    barcodes = '../data/raw_sbm/raw_sbm_data_2024.11.15/barcodes.tsv'
    output_file = '../data/raw_sbm/raw_sbm_data_2024.11.15/adata.h5ad'
    convert_mtx_to_h5ad(mtx, genes, barcodes, output_file)