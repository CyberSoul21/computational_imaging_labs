clear; clc; close all;
% D = load("CI_Lab_NLOS_datasets\Z_d=0.5_l=[1x1]_s=[256x256].mat"); norm_check = 0; confoc_check = 0;
% D = load("CI_Lab_NLOS_datasets\usaf_d=0.5_l=[1x1]_s=[256x256].mat"); norm_check = 1; confoc_check = 0;
% D = load("CI_Lab_NLOS_datasets\bunny_d=0.5_l=[1x1]_s=[256x256].mat"); norm_check = 0; confoc_check = 0;
D = load("CI_Lab_NLOS_datasets\bunny_d=0.5_c=[256x256].mat"); norm_check = 1; confoc_check = 1;

% Resolution
N = 32;

fn = fieldnames(D);
dataset = D.(fn{1});

% % First to-do
% H = squeeze(dataset.data);
% [nx,ny,nt] = size(H);
% 
% x_mid = round(nx/2);
% y_mid = round(ny/2);
% Hx = squeeze(H(:,y_mid,:));
% Hy = squeeze(H(x_mid,:,:));
% 
% figure
% imagesc(Hx)
% axis xy
% colormap hot
% colorbar
% title('x-t slice')
% 
% figure
% imagesc(Hy)
% axis xy
% colormap hot
% colorbar
% title('y-t slice')
% % End first to-do


laserPos = dataset.laserPositions;
spadPos  = dataset.spadPositions;
laserOrigin = dataset.laserOrigin(:)';
spadOrigin  = dataset.spadOrigin(:)';
deltaT = dataset.deltaT;
t0 = dataset.t0;
volCenter = dataset.volumePosition(:)';
volSize = dataset.volumeSize;
H = dataset.data;

if isscalar(volSize)
    volSize = [volSize volSize volSize];
end


% Preprocess spad to avoid loops
disp('preprocess spad');
spadPos_flat = reshape(spadPos, [], 3);
if confoc_check == 0 % No confocal
    [lu,lv,su,sv,nt] = size(H);   
    H_reshaped = reshape(H(1,1,:,:,:), su, sv, nt);
    xl = squeeze(laserPos)';
else % Confocal
    [su,sv,nt] = size(H); %su, sv would be cu, cv in this case   
    H_reshaped = H;      
end

% Grid of voxels
x = linspace(volCenter(1)-volSize(1)/2, volCenter(1)+volSize(1)/2, N);
y = linspace(volCenter(2)-volSize(2)/2, volCenter(2)+volSize(2)/2, N);
z = linspace(volCenter(3)-volSize(3)/2, volCenter(3)+volSize(3)/2, N);
[X,Y,Z] = ndgrid(x,y,z);
G = zeros(N,N,N);

if norm_check == 1 % Normalized
    d1 = 0;
    d4 = 0;
else % Not Normalized
    d1 = norm(laserOrigin - xl);
    d4 = sqrt(sum((spadPos_flat - spadOrigin).^2, 2));  
end

disp('ini reprojection');
for ix = 1:N
    for iy = 1:N
        for iz = 1:N
            xv = [X(ix,iy,iz), Y(ix,iy,iz), Z(ix,iy,iz)];
            
            if confoc_check == 0 % No confocal
                d2 = norm(xl - xv);
                diff = spadPos_flat - xv;                 
                d3 = sqrt(sum(diff.^2, 2));               
                tv = d1 + d2 + d3 + d4;
            else % Confocal
                diff = spadPos_flat - xv;
                d = sqrt(sum(diff.^2, 2));
                tv = 2*d;
            end

            tidx = ceil((tv - t0)/deltaT);
            valid = (tidx > 0) & (tidx <= nt);

            if any(valid)
                tidx_valid = tidx(valid);
                idx_spad = find(valid);

                [iu, iv] = ind2sub([su, sv], find(valid));
                ind = sub2ind([su, sv, nt], iu, iv, tidx_valid);
                G(ix,iy,iz) = sum(H_reshaped(ind));
            else
                G(ix,iy,iz) = 0;
                disp("no valid");
            end
        end
    end
end    
disp('end reprojection');

% This accumulates over all depth
proj = max(G,[],2); %Collapse over plane X-Z
proj = squeeze(proj);
figure
imagesc(rot90(proj, -1))
axis xy
axis image
colormap hot
colorbar
title('Backprojection reconstruction')

%3D representation
volshow(G)
title('Volumen 3D')

% Laplacian reconstruction
f_lap = fspecial3('lap');
G_lap = imfilter(G, -f_lap, 'symmetric');
G_lap(G_lap < 0) = 0;
proj_lap = max(G_lap,[],2); %Collapse over plane X-Z
proj_lap = squeeze(proj_lap);
figure
imagesc(rot90(proj_lap, -1))
axis xy
axis image
colormap hot
colorbar
title('Backprojection + Laplacian')

%laplacian Gauss reconstruction
sigma_v = [0.5, 1.0, 2.0];   % try 0.5, 1, 2
for sigma = sigma_v
    f_log = fspecial3('log', [5 5 5], sigma);
    G_log = imfilter(G, -f_log, 'symmetric');
    G_log(G_log < 0) = 0;
    proj_log = max(G_log,[],2); %Collapse over plane X-Z
    proj_log = squeeze(proj_log); 
    figure
    imagesc(rot90(proj_log, -1))
    axis xy
    axis image
    colormap hot
    colorbar
    title(['Backprojection + LoG (sigma = ', num2str(sigma), ')'])
end