# Assignment 2 - README

## Group:
- Javier Almario 962449 
- Alvaro Provencio 960625

## File Information

- The folder contains files .m, images (results) and the final report as **`CI_Assignment_2.pdf`**.
- First Download the folder: https://drive.google.com/drive/folders/1_dZ6Zsafps7qR4Mcx7pLXrK3ZliQif3R?usp=sharing
- The main files are **`deblurring_demo.m`** and **`power_spect_demo.m`**, which contains the code for Assignment 2.

## How to Use

### 1. **Separate Runs**: 

Each section of the code should be run independently based on which analysis you want to perform. Depending on which part of the code you wish to test (e.g., power spectrum or deblurring), you may need to **comment or uncomment** the relevant lines of code.

To select different apertures (e.g., different window sizes or sigma values), follow the instructions below.

### 2. **Testing with Different Parameters**:
Both scripts are set up to run with **two nested for loops**. These loops allow you to experiment with different combinations of the `sigma` and `blur_size` parameters.

### 2. **Assignment structure**:
project_root/
    ├── power_spect_demo.m
    ├── deblurring_demo.m
    └── results/
        ├── result_sigma_0.5_blur3.png
        ├── result_sigma_1.0_blur5.png
        └── ... (other result files)


### Example:
```matlab
% Nested loop for sigma and blur_size
for sigma = [0.5, 1.0, 1.5, 2.0]  % Example sigma values
    for blur_size = [3, 5, 7]    % Example blur sizes
        % Call the function or processing code here with the current sigma and blur_size
        result = perform_deblurring(input_image, sigma, blur_size);
        display_result(result);  % Display or save result
    end
end

