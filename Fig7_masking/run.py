# # Reproduce Figure 7

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import binary_dilation, binary_erosion, binary_closing
from fastmri.data import transforms, subsample
import torch
import cfl
from bart import bart
import os
import sys
import subprocess

if 'BART_TOOLBOX_PATH' in os.environ and os.path.exists(os.environ['BART_TOOLBOX_PATH']):
	sys.path.append(os.path.join(os.environ['BART_TOOLBOX_PATH'], 'python'))
elif 'TOOLBOX_PATH' in os.environ and os.path.exists(os.environ['TOOLBOX_PATH']):
	sys.path.append(os.path.join(os.environ['TOOLBOX_PATH'], 'python'))
else:
	raise RuntimeError("BART_TOOLBOX_PATH is not set correctly!")


# Utils
# Mask generation function
def create_brain_masks(mri_scan, loose_padding=5, hole_structre=5):
    normalized_scan = (mri_scan - mri_scan.min()) / (mri_scan.max() - mri_scan.min()) # Normalize

    # thresholding
    tight_mask = normalized_scan > 0.1
    # Fill holes in the tight mask using binary closing
    structure = np.ones((hole_structre, hole_structre))
    structure_erosion = np.ones((10, 10))

    tight_mask = binary_closing(tight_mask, structure=structure)
    tight_mask = binary_erosion(tight_mask, structure=structure_erosion)

    # Create a wide mask by applying dilation to the tight mask
    loose_mask = binary_dilation(tight_mask, structure=structure, iterations=loose_padding)
    loose_mask = binary_erosion(loose_mask, structure=structure_erosion)

    return tight_mask.astype(np.uint8), loose_mask.astype(np.uint8)

def ifft2c(kspace):
    return np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(kspace, axes=(0, 1)), axes=(0, 1)), axes=(0, 1))

def fft2c(image):
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(image, axes=(0, 1)), axes=(0, 1)), axes=(0, 1))

def calc_NRMSE(I_pred,I_true):
        # Reshape the images into vectors
        I_true = np.reshape(I_true,(1,-1))
        I_pred = np.reshape(I_pred,(1,-1))
        # Mean Square Error
        MSE = np.square(np.subtract(I_true,I_pred)).mean()
        # Root Mean Square Error
        RMSE = np.sqrt(MSE)
        # Normalized Root Mean Square Error
        rr = np.max(I_true) - np.min(I_true) # range
        NRMSE = RMSE/rr
        return NRMSE

def crop_image_vertical(data, crop_fraction=0.2):
    """
    Crops the upper and lower parts of the data.
    Args:
        data: 2D or 3D array to crop.
        crop_fraction: Fraction to crop from the top and bottom.
    Returns:
        Cropped data.
    """
    height = data.shape[0]
    crop_size = int(height * crop_fraction)
    return data[crop_size:height - crop_size, ...]

def crop_image_horizontal(data, crop_fraction=0.2):
    """
    Crops the upper and lower parts of the data.
    Args:
        data: 2D or 3D array to crop.
        crop_fraction: Fraction to crop from the top and bottom.
    Returns:
        Cropped data.
    """
    width = data.shape[1]
    crop_size = int(width * crop_fraction)
    return data[:,crop_size:width - crop_size, ...]

# %% [markdown]
# # Read raw data

# %%

# %%
# local data
odir = os.path.abspath(os.path.dirname(sys.argv[0]))
ksp_all = cfl.readcfl(odir+'/../Data/ksp_fully')
print(ksp_all.shape)
slice_ksp_r1 = ksp_all
print(slice_ksp_r1.shape)
slice_ksp_r1 = slice_ksp_r1[:,:,0,:]
print(slice_ksp_r1.shape)

# %% [markdown]
# # Generate sampling mask

# %%
def get_mask_func( factor):
    center_fractions = 0.06 * 4/factor # EquiSpacedMaskFunc
    mask_func = subsample.EquiSpacedMaskFunc(
    center_fractions=[center_fractions],
    accelerations=[factor],
    )
    return mask_func
mask_func = get_mask_func(2)
mask = transforms.apply_mask(torch.ones(640,420,32), mask_func)[0]
undersampling_mask = mask.numpy()

# %% [markdown]
# # Data processing

# %%
# Cuts of background
slice_vertical = 0.26
slice_horizontal = 0.2

#Prewhiten
scan = bart(1,'fft -i -u 1',slice_ksp_r1 )
noise = bart(1,'transpose 1 0',scan)
ksp_white = bart(1,'whiten -n',slice_ksp_r1, noise)

slice_ksp_r1 = ksp_white  # the data is kspace
slice_r1 = ifft2c(slice_ksp_r1)


slice_rss_r1 = np.sqrt(np.sum(np.abs(slice_r1) ** 2, axis=-1)) # RSS of scan
slice_rss_r1 = np.flipud(slice_rss_r1[:,:]) # Flipup
slice_rss_r1_cropped = crop_image_vertical(slice_rss_r1, crop_fraction=slice_vertical)  # Crop the RSS
slice_rss_r1_cropped = crop_image_horizontal(slice_rss_r1_cropped, crop_fraction=slice_horizontal)  # Crop the RSS

# Generate accelerated data using the fully-sampled data
slice_ksp_r2 = slice_ksp_r1 * undersampling_mask
slice_r2 = ifft2c(slice_ksp_r2)
slice_rss_r2 = np.sqrt(np.sum(np.abs(slice_r2) ** 2, axis=-1)) # RSS of scan
slice_rss_r2 = np.flipud(slice_rss_r2[:,:]) # Flipup
slice_rss_r2_cropped = crop_image_vertical(slice_rss_r2, crop_fraction=slice_vertical)  # Crop the RSS
slice_rss_r2_cropped = crop_image_horizontal(slice_rss_r2_cropped, crop_fraction=slice_horizontal)  # Crop the RSS

r2_image = slice_rss_r2_cropped
r1_image = slice_rss_r1_cropped

# %% [markdown]
# # Metric scores calculation masks generation

# %%
#Masks generation
tight_mask, wide_mask = create_brain_masks(r1_image, loose_padding=1,hole_structre=40)
no_mask = np.ones(r1_image.shape) # Only frame
no_mask[:,-3:-1] = 0
no_mask[:,0:3] = 0
no_mask[0:3,:] = 0
no_mask[-3:-1,:] = 0


# %% [markdown]
# # Calculate NRMSE scores

# %%
# NRMSE calculation
no_mask_nrmse = calc_NRMSE(r2_image,r1_image)
wide_mask_nrmse = calc_NRMSE(r2_image*wide_mask,r1_image*wide_mask)
tight_mask_nrmse = calc_NRMSE(r2_image*tight_mask,r1_image*tight_mask)

print(f'NRMSE with no mask: {no_mask_nrmse:.3f}')
print(f'NRMSE with loose mask: {wide_mask_nrmse:.3f}')
print(f'NRMSE with tight mask: {tight_mask_nrmse:.3f}')


# # Combine results with inkscape template for final figure

## For Saving plot as the paper

# Create plots of Fully-Sampled & Reconstruction
plt.figure(figsize=(5, 8))
plt.imshow(np.sqrt(r2_image), cmap='gray')
plt.axis('off')
# Save as PNG
plt.savefig("plot_2.png", dpi=900, bbox_inches='tight', pad_inches=0)
plt.close()

plt.figure(figsize=(5, 5))
plt.imshow(np.sqrt(r1_image), cmap='gray')
plt.axis('off')
plt.savefig("plot_1.png", dpi=900, bbox_inches='tight', pad_inches=0)
plt.close()

# Masks plots
# Add text in the middle of the plot
middle_x = (np.sqrt(no_mask)).shape[1] // 2  # Middle x-coordinate
middle_y = (np.sqrt(no_mask)).shape[0] // 2  # Middle y-coordinate 

plt.figure(figsize=(5, 5))
no_mask_label = f"NRMSE: {no_mask_nrmse:.3f}".replace("0.", ".")
plt.title(no_mask_label, fontsize=32)
plt.imshow(np.sqrt(no_mask), cmap='gray')
plt.text(middle_x, middle_y, "No \nmask", fontsize=32, color='black', ha='center', va='center')
plt.axis('off')
plt.savefig("plot_3.png", dpi=900, bbox_inches='tight', pad_inches=0)
plt.close()

plt.figure(figsize=(5, 5))
wide_mask_label = f"NRMSE: {wide_mask_nrmse:.3f}".replace("0.", ".")
plt.title(wide_mask_label, fontsize=32)
plt.imshow(np.sqrt(wide_mask), cmap='gray')
plt.text(middle_x, middle_y, "Wide \nmask", fontsize=32, color='black', ha='center', va='center')
plt.axis('off')
plt.savefig("plot_4.png", dpi=900, bbox_inches='tight', pad_inches=0)
plt.close()

plt.figure(figsize=(5, 5))
tight_mask_label = f"NRMSE: {tight_mask_nrmse:.3f}".replace("0.", ".")
plt.title(tight_mask_label, fontsize=32)
plt.imshow(np.sqrt(tight_mask), cmap='gray')
plt.text(middle_x, middle_y, "Tight \nmask", fontsize=32, color='black', ha='center', va='center')
plt.axis('off')
plt.savefig("plot_5.png", dpi=900, bbox_inches='tight', pad_inches=0)
plt.close()


subprocess.run(["inkscape", "--export-filename=Fig_7_NRMSE_vs_mask.pdf", "Fig_7_template.svg"])
#


