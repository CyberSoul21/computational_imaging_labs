# Assignment 1 - README

## Group:
- Javier Almario 962449 
- Alvaro Provencio 960625


This project demonstrates a sequence of image processing steps using MATLAB. The main focus is on reading and processing images, including linearization, demosaicing, white balancing, denoising, color balancing, tone reproduction, and image compression.

## File Information

- The folder contains files .mlx, images (results) and the final report as **`CI_Assignment_1.pdf`**.
- First Download the folder: https://drive.google.com/drive/folders/1Dwoza37cJ_yw2oLVDJAhLAj7Thr-KH0k?usp=sharing
- The main file to execute is **`lab1.mlx`**, which contains the requirements for Assignment 1.
- **`lab1_Extra.mlx`** is a complete pipeline with additional image and processing steps.

## Folder Structure

- **input_data**: Contains the input image files required for the processing pipeline.
    - `cornell_box.tiff`: The TIFF image used for the workflow.
    - `cornell_box.CR2`: A raw image file used for demosaicing.
- **output_lab1**: Contains the processed images saved during various stages of the workflow.
    - `02_data.png`: Processed version of the input image with intensity scaling.
    - `02_double.png`: Image after reading and converting to double precision.
    - `02_linearized.png`: Linearized image after dark noise and saturation adjustments.
    - `03_NN_interpolation.png`: Demosaiced image using nearest neighbor interpolation.
    - `03_Bilinear_RGGB.png`: Demosaiced image using bilinear interpolation.
    - `04_white_manual.png`: Image after manual white balance.
    - `04_white_world.png`: Image after automatic white world white balance.
    - `04_gray_world.png`: Image after automatic gray world white balance.
    - `05_mean_denoising.png`: Image after applying mean filtering denoising.
    - `05_median_denoising.png`: Image after applying median filtering denoising.
    - `05_gaussian_denoising.png`: Image after applying Gaussian filtering denoising.
    - `06_color_balance.png`: Image after applying color balance.
    - `07_final_tone.png`: Final image after tone reproduction adjustments.
    - `08_final_tone.png`: Final image after compression and tone reproduction adjustments.
    - `08_final_tone.jpg`: JPEG version of the final image.

- **output_extra**: Contains additional images produced from the full pipeline in `lab1_Extra.mlx`. The naming convention for the images is similar to those in **`output_lab1`**, with additional processing steps applied.
    - The images in this folder follow the same naming convention but are based on other input images and extended processes described in **`lab1_Extra.mlx`**.

