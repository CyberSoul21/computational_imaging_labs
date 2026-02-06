function f = zDefocused(f0, k, sigma, flag)

	% [height, width] = size(f0);
	[height, width, channel] = size(f0); %Exercise 5
    f = zeros([height, width, channel]);

    if flag == 1
        % Assume the defocus is circular convolution
		k = zPSFPad(k, height, width);
		k = fft2(fftshift(k));
        % f = abs(ifft2(fft2(f0) .* k)) + randn(height, width) * sigma;
        %Exercise 5
        for c = 1:channel
		    f(:,:,c) = abs(ifft2(fft2(f0(:,:,c)) .* k)) + randn(height, width) * sigma;
        end
    else
        % Not make that assumption. But in this case, because of the 
        % boundary effects, the "effective" noise level is very high.
		k = rot90(k, 2);
		% f = imfilter(f0, k) + randn(height, width) * sigma;
        %Exercise 5
        for c = 1:channel
		    f(:,:,c) = imfilter(f0(:,:,c), k) + randn(height, width) * sigma;
        end
        
    end

end
