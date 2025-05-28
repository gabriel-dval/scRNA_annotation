'''
Utility file to make some helpful functions across all annotation
methods.
'''

import anndata
import pandas as pd
import numpy as np
from PIL import Image
import cv2

import os
from scipy.io import mmread, mmwrite

# Utility functions

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

    return adata


def split_csv_dataset(path_to_file, number_of_parts, output_dir):
    '''Function to split a csv file into multiple parts and save the
    resulting files in the same directory

    Args
    ----
    path_to_file : str
        The path to the csv file to split
    number_of_parts : int
        The number of parts to split the file into
    output_dir : str
        The directory to save the resulting files

    Returns
    -------
    None
    '''
    # Read the file
    data = pd.read_csv(path_to_file)
    data.set_index('Unnamed: 0', inplace=True)

    # Create the output directory
    os.makedirs(output_dir, exist_ok=True)

    # Calculate the number of columns per part
    columns_per_part = len(data.columns) // number_of_parts
    for i in range(number_of_parts):
        start_col = i * columns_per_part
        end_col = (i + 1) * columns_per_part if i != number_of_parts - 1 else len(data.columns)
        part = data.iloc[:, start_col:end_col]
        print(part.head())
        part.to_csv(os.path.join(output_dir, f'log_immune_norm_batch{i + 1}.csv'), index=True)


def convert_nonblack_to_white(image_path, output_path, threshold=50):
    # Load the image
    image = Image.open(image_path)
    
    # Convert image to RGB mode
    image = image.convert("RGB")
    
    # Get pixel data
    pixels = image.load()
    width, height = image.size
    
    # Process each pixel
    # Process each pixel
    for x in range(width):
        for y in range(height):
            r, g, b = pixels[x, y]
            # Calculate brightness (approximate luminance)
            brightness = (r + g + b) / 3
            # If the pixel is not dark enough, turn it white
            if brightness > threshold:
                pixels[x, y] = (255, 255, 255)
    
    # Save the processed image
    image.save(output_path)


def convert_nonblack_to_transparent(image_path, output_path, threshold=50):
    # Load the image
    image = Image.open(image_path)
    
    # Convert image to RGBA mode (adds alpha channel)
    image = image.convert("RGBA")
    
    # Get pixel data
    pixels = image.load()
    width, height = image.size
    
    # Process each pixel
    for x in range(width):
        for y in range(height):
            r, g, b, a = pixels[x, y]
            # Calculate brightness (approximate luminance)
            brightness = (r + g + b) / 3
            # If the pixel is bright (white or near white), make it transparent
            if brightness > threshold:
                pixels[x, y] = (255, 255, 255, 0)  # Set alpha to 0 (transparent)
            # If pixel is dark, set it to white
            if brightness < threshold:
                pixels[x, y] = (255, 255, 255, a)  # Set to white, preserve alpha
            
    
    # Save the processed image
    image.save(output_path, "PNG")  # Save as PNG to preserve transparency


def thicken_lines(input_path, output_path, thickness=3):
    # Load image with alpha channel
    img = Image.open(input_path).convert("RGBA")
    img_np = np.array(img)
    
    # Separate channels
    r, g, b, a = img_np[:, :, 0], img_np[:, :, 1], img_np[:, :, 2], img_np[:, :, 3]
    
    # Create a mask where white lines are (RGB = 255 and alpha > 0)
    white_mask = (r == 255) & (g == 255) & (b == 255) & (a > 0)
    
    # Convert mask to uint8 image for OpenCV
    mask_uint8 = np.uint8(white_mask) * 255

    # Dilate (thicken) the white lines
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (thickness * 2 + 1, thickness * 2 + 1))
    dilated_mask = cv2.dilate(mask_uint8, kernel)

    # Create a new image with transparent background
    new_img_np = np.zeros_like(img_np)
    
    # Set white color where dilated mask is present
    new_img_np[dilated_mask > 0] = [255, 255, 255, 255]  # white with full opacity
    
    # Convert back to image and save
    new_img = Image.fromarray(new_img_np, mode="RGBA")
    new_img.save(output_path)


def white_to_transparent(image_path, output_path=None, threshold=240):
    """
    Convert white or near-white pixels to transparent.
    
    Args:
        image_path (str): Path to input image
        output_path (str, optional): Path for output image. If None, adds '_transparent' to original name
        threshold (int): RGB threshold value (0-255). Pixels with all RGB values >= threshold become transparent
    
    Returns:
        PIL.Image: Image with transparent background
    """
    # Open image and convert to RGBA if not already
    img = Image.open(image_path).convert("RGBA")
    
    # Convert to numpy array for easier manipulation
    data = np.array(img)
    
    # Find pixels where all RGB values are >= threshold (near-white)
    white_pixels = np.all(data[:, :, :3] >= threshold, axis=2)
    
    # Set alpha channel to 0 (transparent) for white pixels
    data[white_pixels, 3] = 0
    
    # Convert back to PIL Image
    result = Image.fromarray(data)
    
    # Save if output path provided
    if output_path:
        result.save(output_path)
    elif output_path is None and isinstance(image_path, str):
        # Auto-generate output filename
        name, ext = image_path.rsplit('.', 1)
        result.save(f"{name}_transparent.png")
    
    return result



if __name__ == '__main__':
    mtx = '../data/raw_sbm/raw_sbm_data_2024.11.15/matrix.mtx'
    genes = '../data/raw_sbm/raw_sbm_data_2024.11.15/genes.tsv'
    barcodes = '../data/raw_sbm/raw_sbm_data_2024.11.15/barcodes.tsv'
    output_file = '../data/raw_sbm/raw_sbm_data_2024.11.15/adata.h5ad'
    #convert_mtx_to_h5ad(mtx, genes, barcodes, output_file)

    #Test extract_gem_data
    ref_file = '../../Documents/data/ref_sbm/BoneMarrow_cellxgene.h5ad'
    #res = extract_gem_data(ref_file, '../data/test_datasets/mATLAS_Brain_Myeloid_facs')
    adata = anndata.read_h5ad(ref_file)
    print(adata)
    #Filter the AnnData object to keep only cells from 'Bone Marrow'
    #Now extract only this bone marrow data
    # res = extract_gem_data(ref_file, '../data/test_datasets/Azimuth_HLA')
    # adata.obs['ann_level_3'].to_csv(os.path.join('../data/test_datasets/Azimuth_HLA', 
    #                                          "cluster_names.tsv"), 
    #                                          sep="\t", index=False, header=False)
    # adata.obs['leiden_3'].to_csv(os.path.join('../data/test_datasets/Azimuth_HLA', 
    #                                          "cluster_numbers.tsv"), 
    #                                          sep="\t", index=False, header=False)
    

    # Split the dataset
    # split_csv_dataset('../data/raw_sbm/log_immune_norm.csv', 
    #                   4, 
    #                   '../data/raw_sbm/split_data')
    
    #data = pd.read_csv( '../data/raw_sbm/split_data/log_immune_norm_batch1.csv')
    #print(data.loc[:,'Unnamed: 0'])

    # for i in range(1):
    #     convert_nonblack_to_white(f'../../Desktop/{i+1}.jpeg', f'../../Desktop/{i+1}_white.jpeg')
    #     convert_nonblack_to_transparent(f'../../Desktop/{i+1}_white.jpeg', f'../../Desktop/{i+1}_transparent.png', )
    #     thicken_lines(f'../../Desktop/{i+1}_transparent.png', f'../../Desktop/{i+1}_transparent_thickened.png', 3)
 
    white_to_transparent('../../Desktop/g14.png','../../Desktop/mc_face_transparent.png')
    