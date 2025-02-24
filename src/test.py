'''
Utility file to make some helpful functions across all annotation
methods.
'''

import anndata
import pandas as pd

import os
from scipy.io import mmread, mmwrite

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


def extract_gem_data(adata, output_dir=None):
    '''
    Function to extract the GEM data from the adata file and corresponding
    barcode and genes names.

    Parameters
    ---
    adata : anndata.AnnData
        The AnnData object to extract the data from.
    output_dir : str
        The output directory to save the extracted data.

    Returns
    ---
    mtx :
        Sparse matrix of the data (able to be loaded in R)
    genes :
        List of gene names
    barcodes :
        List of cell barcodes
    '''
    # First read adata using anndata
    adata = anndata.read_h5ad(adata)

    # Create results directory
    os.makedirs(output_dir, exist_ok=True)

    # Extract the matrix, genes and barcodes
    mtx = adata.X.T # Transpose the matrix for 10x format
    mmwrite(os.path.join(output_dir, "matrix.mtx"), mtx)

    # Genes
    genes = pd.DataFrame(adata.var_names)  # Gene names
    genes.to_csv(os.path.join(output_dir, "genes.tsv"), sep="\t", index=False, header=False)

    # Barcodes
    barcodes = pd.DataFrame(adata.obs_names)  # Cell barcodes
    barcodes.to_csv(os.path.join(output_dir, "barcodes.tsv"), sep="\t", index=False, header=False)

    # Also extract and save annotation labels
    if 'cluster_names' in adata.obs.columns:
        adata.obs['cluster_names'].to_csv(os.path.join(output_dir, "cluster_names.tsv"), 
                                          sep="\t", index=False, header=False)

    return (mtx, genes, barcodes)


if __name__ == '__main__':
    mtx = '../data/raw_sbm/raw_sbm_data_2024.11.15/matrix.mtx'
    genes = '../data/raw_sbm/raw_sbm_data_2024.11.15/genes.tsv'
    barcodes = '../data/raw_sbm/raw_sbm_data_2024.11.15/barcodes.tsv'
    output_file = '../data/raw_sbm/raw_sbm_data_2024.11.15/adata.h5ad'
    #convert_mtx_to_h5ad(mtx, genes, barcodes, output_file)

    # Test extract_gem_data
    ref_file = '../data/ref_sbm/mATLAS_Marrow_droplet.h5ad'
    res = extract_gem_data(ref_file, '../data/ref_sbm//mATLAS_Marrow_droplet_10x')
 
  
    