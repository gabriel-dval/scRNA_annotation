"""
Combined image processor that can:
1. Convert non-black pixels to white or transparent
2. Optionally thicken lines
3. Run from command line with input/output file arguments
"""

import numpy as np
import cv2
from PIL import Image


def convert_nonblack_to_transparent(image_path, output_path, black_to_white = True,threshold=50):
    '''What this does : 
    Converts non-black pixels in an image to transparent, 
    while turning dark pixels to white optionally (black_to_white option).

    Args
    ----
    image_path (str): Path to the input image
    output_path (str): Path to save the processed image
    black_to_white (bool): If True, dark pixels are set to white; if False
    they remain unchanged
    threshold (int): Brightness threshold (0-255) for determining what's "black"
    (default is 50, which is a low threshold to catch most dark pixels)

    Returns
    -------
    None: The processed image is saved to output_path.
    '''
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
            if brightness < threshold and black_to_white:
                pixels[x, y] = (255, 255, 255, a)  # Set to white, preserve alpha
            
    
    # Save the processed image
    image.save(output_path, "PNG")  # Save as PNG to preserve transparency


def thicken_lines(input_path, output_path, thickness=3):
    """
    Thicken white lines in an image by applying morphological dilation.

    Args
    ----
        input_path (str): Path to input image
        output_path (str): Path to save processed image
        thickness (int): Thickness factor for line dilation (default is 3)
    
    Returns
    -------
        None: The processed image is saved to output_path.
    """
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
    
    Args
    ----
        image_path (str): Path to input image
        output_path (str, optional): Path for output image. If None, adds '_transparent' to original name
        threshold (int): RGB threshold value (0-255). Pixels with all RGB values >= threshold become transparent
    
    Returns
    -------
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
    
    # Example
    convert_nonblack_to_transparent('image.png','image_changed.png')