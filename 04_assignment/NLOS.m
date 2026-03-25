clear; clc; close all;
D = load("CI_Lab_NLOS_datasets\Z_d=0.5_l=[1x1]_s=[256x256].mat"); norm_check = 1;
% D = load("CI_Lab_NLOS_datasets\usaf_d=0.5_l=[1x1]_s=[256x256].mat"); norm_check = 0;
% D = load("CI_Lab_NLOS_datasets\bunny_d=0.5_l=[1x1]_s=[256x256].mat"); norm_check = 1;

fn = fieldnames(D);
dataset = D.(fn{1});

% extraer histograma
H = squeeze(dataset.data);

% First to-do
[nx,ny,nt] = size(H);

% slices centrales
x_mid = round(nx/2);
y_mid = round(ny/2);

Hx = squeeze(H(:,y_mid,:));
Hy = squeeze(H(x_mid,:,:));

% mejorar visibilidad
Hx = log(Hx + 1);
Hy = log(Hy + 1);

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

% End first to-do


laserPos = dataset.laserPositions;
spadPos  = dataset.spadPositions;
% laserOrigin = dataset.laserOrigin;
% spadOrigin  = dataset.spadOrigin;
laserOrigin = dataset.laserOrigin(:)';
spadOrigin  = dataset.spadOrigin(:)';
deltaT = dataset.deltaT;
t0 = dataset.t0;
% volCenter = dataset.volumePosition;
volCenter = dataset.volumePosition(:)';
volSize = dataset.volumeSize;

if isscalar(volSize)
    volSize = [volSize volSize volSize];
end

% Resolution
N = 32;

% Histograma (NO squeeze)
H = dataset.data;
[lu,lv,su,sv,nt] = size(H);


% Grid of voxels
x = linspace(volCenter(1)-volSize(1)/2, volCenter(1)+volSize(1)/2, N);
y = linspace(volCenter(2)-volSize(2)/2, volCenter(2)+volSize(2)/2, N);
z = linspace(volCenter(3)-volSize(3)/2, volCenter(3)+volSize(3)/2, N);
[X,Y,Z] = meshgrid(x,y,z);
G = zeros(N,N,N);

% Preprocess spad to avoid loops
disp('preprocess spad');
Ns = su * sv;
spadPos_flat = reshape(spadPos, [], 3);   
H_reshaped = reshape(H(1,1,:,:,:), su, sv, nt);

xl = squeeze(laserPos(1,1,:))';
if norm_check == 1
    d1 = norm(laserOrigin - xl);
    d4 = sqrt(sum((spadPos_flat - spadOrigin).^2, 2));  
else
    d1 = 0;
    d4 = 0;
end

disp('ini reprojection');
for ix = 1:N
    for iy = 1:N
        for iz = 1:N
            xv = [X(ix,iy,iz), Y(ix,iy,iz), Z(ix,iy,iz)];
            d2 = norm(xl - xv);

            diff = spadPos_flat - xv;                 
            d3 = sqrt(sum(diff.^2, 2));               

            tv = d1 + d2 + d3 + d4;
            tidx = ceil((tv - t0)/deltaT);
        
            valid = (tidx > 0) & (tidx <= nt);

            if any(valid)
                tidx_valid = tidx(valid);
                idx_spad = find(valid);
        
                [iu, iv] = ind2sub([su, sv], idx_spad);
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
proj = max(G,[],1); %Collapse over plane Y-Z
proj = squeeze(proj);
% proj = log(proj + 1); %Improve visualization
figure
% imagesc(proj)
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
proj_lap = max(G_lap,[],1);
proj_lap = squeeze(proj_lap);
figure
% imagesc(proj_lap)
imagesc(rot90(proj_lap, -1))
% imagesc(log(proj_lap+1))
axis xy
axis image
colormap hot
colorbar
title('Backprojection + Laplacian')

%laplacian Gauss reconstruction
sigma = 1.0;   % try 0.5, 1, 2
f_log = fspecial3('log', [5 5 5], sigma);
G_log = imfilter(G, -f_log, 'symmetric');
G_log(G_log < 0) = 0;
proj_log = max(G_log,[],1);
proj_log = squeeze(proj_log);
figure
% imagesc(proj_log)
imagesc(rot90(proj_log, -1))
% imagesc(log(proj_log+1))
axis xy
axis image
colormap hot
colorbar
title(['Backprojection + LoG (sigma = ', num2str(sigma), ')'])