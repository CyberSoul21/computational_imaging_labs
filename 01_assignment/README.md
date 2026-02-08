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
    - `cornell_box_2.png`: Linearized image after dark noise and saturation adjustments.
    - `cornell_box_3.png`: Demosaiced image using nearest neighbor interpolation RGGB.
    - `cornell_box_4.png`: Demosaiced image using nearest neighbor interpolation BGGR.
    - `cornell_box_5.png`: Demosaiced image using bilinear interpolation RGGB.
    - `cornell_box_6.png`: Image after automatic gray world white balance.
    - `cornell_box_7.png`: Image after automatic white world white balance.
    - `cornell_box_8.png`: Image after manual white balance.
    - `cornell_box_9.png`: Image before applying denoising.
    - `cornell_box_10.png`: Image after applying Gaussian filtering denoising.
    - `cornell_box_11.png`: Image after applying mean filtering denoising.
    - `cornell_box_12.png`: Image after applying median filtering denoising.
    - `cornell_box_13.png`: Image after applying color balance.
    - `cornell_box_14.png`: Final image after tone reproduction adjustments.
    - `cornell_box_15.jpg`: JPEG version of the final image with 5.28 compression ratio.
    - `cornell_box_16.png`: PNG version of the final image.
    - `cornell_box_17.jpg`: JPEG version of the final image with 58.03 compression ratio.

- **output_extra**: Contains additional images produced from the full pipeline in `lab1_Extra.mlx`. The naming convention for the images is similar to those in **`output_lab1`**, with additional processing steps applied.
    - The images in this folder follow the same naming convention but are based on other input images and extended processes described in **`lab1_Extra.mlx`**.

