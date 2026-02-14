function outK = zPSFPad(inK, height, width)
    % This function is to zeropadding the psf
    
	[sheight, swidth] = size(inK);

	outK = zeros(height, width);
	
    outK(floor(end/2-sheight/2) + 1:floor(end/2-sheight/2) + sheight, ...
         floor(end/2-swidth/2) + 1:floor(end/2-swidth/2) + swidth) = inK;

end

% Read data
% aperture = imread('apertures/circular.bmp');
% aperture = imread('apertures/Levin.bmp');
% aperture = imread('apertures/raskar.bmp');
aperture = imread('apertures/zhou.bmp');
image = imread('images/penguins.jpg');

% Exercise 5
image = image(:, :, 1); %Comment for color image, Uncomment for gray

% Noise level (Gaussian noise)
% sigma_v = 0.005;
%sigma_v = [0.0005, 0.005, 0.010, 0.1, 1, 2];
sigma_v = [0.001, 0.01, 0.05];

% Blur size
%blurSize = 7;
%blurSize_v = [1,7,14];
blurSize_v = [5, 13, 25];



f0 = im2double(image);
[height, width, channel] = size(f0);

for sigma = sigma_v
    for blurSize = blurSize_v
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
        
        % Recover
        psf = zPSFPad(k1, height, width);
        f0_hat = zDeconvWNR(f1, k1, C); % Wienner with priors
        % f0_hat = deconvlucy(f1, psf); %Lucy
        % f0_hat = deconvwnr(f1, psf); %Wienner without priors
        
        % Display results
        fig = figure;
        
        subplot_tight(1, 3, 1, 0.0, false)
        imshow(f0);
        title('Focused');
        
        subplot_tight(1, 3, 2, 0.0, false)
        imshow(f1);
        title('Defocused');
        
        subplot_tight(1, 3, 3, 0.0, false)
        imshow(f0_hat);
        title_s = strcat('Recovered' , ' Sigma: ' , string(sigma), ' BlurSize: ',string(blurSize));
        title(title_s);

        save_name = strcat('results/exercise2/zhou/sigma_',string(sigma),'_blursize_',string(blurSize),'.png');
        saveas(fig,save_name)
    end
end




