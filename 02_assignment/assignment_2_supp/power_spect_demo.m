% Read data
aperture = imread('apertures/circular.bmp');
image = imread('images/penguins.jpg');
image = image(:, :, 1);

% Noise level (gaussian noise)
sigma = 0.005;

% Blur size
blurSize = 7;


f0 = im2double(image);
[height, width, channel] = size(f0);

% Prior matrix: 1/f law
A_star = eMakePrior(height, width) + 0.00000001;
C = sigma.^2 * height * width ./ A_star;

% Normalization
temp = fspecial('disk', blurSize);
flow = max(temp(:));

% Calculate effective PSF
k1 = im2double(...
    imresize(aperture, [2*blurSize + 1, 2*blurSize + 1], 'nearest')...
);

k1 = k1 * (flow / max(k1(:)));

% Apply blur
f1 = zDefocused(f0, k1, sigma, 0);

% Padding aperture
k1P = zPSFPad(k1, max(height, width), max(height, width));

% Aperture power spectra
F = fft2(k1P);
F = fftshift(F .* conj(F));

S = log(F);
S_X = S(:, round(length(S)/2) + 1);
S_Y = S(round(length(S)/2) + 1, :);

% Display results
figure;

subplot_tight(2, 2, 1, 0.05, false)
imagesc(k1P);
axis('image')
axis('off')
title('Aperture');

subplot_tight(2, 2, 2, 0.05, false)
imagesc(S);
axis('image')
axis('off')
title('Aperture Frequency');

subplot_tight(2, 2, 3, 0.05, false)
plot(linspace(-1, 1, length(S_X)), S_X)
grid('on')
title('Normalized frecuency X');

subplot_tight(2, 2, 4, 0.05, false)
plot(linspace(-1, 1, length(S_Y)), S_Y)
grid('on')
title('Normalized frecuency Y');

% Image power spectra
F0 = fft2(f0);
F0 = fftshift(F0 .* conj(F0));

F1 = fft2(f1);
F1 = fftshift(F1 .* conj(F1));

S0 = log(F0);
S0_X = S0(:, round(length(S0)/2) + 1);
S0_Y = S0(round(length(S0)/2) + 1, :);

S1 = log(F1);
S1_X = S1(:, round(length(S1)/2) + 1);
S1_Y = S1(round(length(S1)/2) + 1, :);

% Display results
figure;

subplot_tight(2, 2, 1, [0.1 0.05], false)
plot(linspace(-1, 1, length(S0_X)), S0_X)
grid('on')
title('Original X');

subplot_tight(2, 2, 2, [0.1 0.05], false)
plot(linspace(-1, 1, length(S0_Y)), S0_Y)
grid('on')
title('Original Y');

subplot_tight(2, 2, 3, [0.1 0.05], false)
plot(linspace(-1, 1, length(S1_X)), S1_X)
grid('on')
title('Defocused X');

subplot_tight(2, 2, 4, [0.1 0.05], false)
plot(linspace(-1, 1, length(S1_Y)), S1_Y)
grid('on')
title('Defocused Y');


function outK = zPSFPad(inK, height, width)
    % This function is to zeropadding the psf

    [sheight, swidth] = size(inK);

    outK = zeros(height, width);

    outK(floor(end/2-sheight/2) + 1:floor(end/2-sheight/2) + sheight, ...
         floor(end/2-swidth/2) + 1:floor(end/2-swidth/2) + swidth) = inK;

end
