clc; clear; close all;

%% 1. Load Dataset
% Ensure the path to your dataset folder is correct
dataFolder = '/home/javier/Documents/Computational Imaging/CI_Lab_NLOS_datasets';
fileName = 'Z_d=0.5_l=[1x1]_s=[256x256].mat';
S = load(fullfile(dataFolder, fileName));
D = S.data;

H = squeeze(double(D.data));            % [Nx, Ny, Nt]
spadPos = double(D.spadPositions);      % [Nx, Ny, 3]
laserPos = double(reshape(D.laserPositions, [1, 1, 3])); 
laserOrigin = double(D.laserOrigin(:));
spadOrigin = double(D.spadOrigin(:));
volPos = double(D.volumePosition(:));
volSize = double(D.volumeSize(:));
deltaT = double(D.deltaT);

[NxWall, NyWall, Nt] = size(H);

%% 2. Setup Voxel Grid (Section 2.2: Resolution)
Nx = 64; Ny = 64; Nz = 64; 
x_v = linspace(volPos(1), volPos(1) + volSize(1), Nx);
y_v = linspace(volPos(2), volPos(2) + volSize(2), Ny);
z_v = linspace(volPos(3), volPos(3) + volSize(3), Nz);
[Xv, Yv, Zv] = meshgrid(x_v, y_v, z_v);

G = zeros(Nx, Ny, Nz);

%% 3. Pre-calculate Fixed Distances
% Distance 1: Laser Device to Laser Wall Point
xl = reshape(laserPos, [1, 3]);
d1 = norm(laserOrigin - xl');

% Distance 2: Laser Wall Point to every Voxel (Vectorized)
d2 = sqrt((Xv - xl(1)).^2 + (Yv - xl(2)).^2 + (Zv - xl(3)).^2);

%% 4. Backprojection (Subsampled for speed)
stepWall = 4; % Use 1/4 of the wall points for better quality than step 8
tic;
for sx = 1:stepWall:NxWall
    for sy = 1:stepWall:NyWall
        % SPAD Wall Point
        xs = [spadPos(sx, sy, 1), spadPos(sx, sy, 2), spadPos(sx, sy, 3)];
        
        % Distance 4: SPAD Wall Point to Detector Device
        d4 = norm(xs - spadOrigin');
        
        % Distance 3: Voxel to SPAD Wall Point (Vectorized over whole volume)
        d3 = sqrt((Xv - xs(1)).^2 + (Yv - xs(2)).^2 + (Zv - xs(3)).^2);
        
        % Total Travel Time (Distance)
        dist_total = d1 + d2 + d3 + d4;
        
        % Convert to Bins
        bins = round(dist_total / deltaT) + 1;
        
        % Accumulate only valid bins
        mask = (bins >= 1 & bins <= Nt);
        G(mask) = G(mask) + H(sx, sy, bins(mask));
    end
    if mod(sx, 20) == 1, fprintf('Processing wall row %d/%d...\n', sx, NxWall); end
end
reconstructionTime = toc;
fprintf('Reconstruction finished in %.2f seconds.\n', reconstructionTime);

%% 5. Filtering (Section 2.2: Filtering the reconstruction)
% Use LoG to handle the noise seen in your previous output
h_log = fspecial3('log', [5 5 5], 0.5); 
G_filtered = imfilter(G, -h_log, 'symmetric');
G_filtered(G_filtered < 0) = 0; 

%% 6. Visualization (Section 2.2: 2D Projections)
mip_unfiltered = max(G, [], 3);
mip_filtered = max(G_filtered, [], 3);

figure('Position', [100, 100, 1000, 400]);
subplot(1,2,1);
imagesc(x_v, y_v, mip_unfiltered'); axis image; colormap hot; colorbar;
title(['Unfiltered MIP (Time: ', num26str(reconstructionTime, '%.1f'), 's)']);
xlabel('X (m)'); ylabel('Y (m)');

subplot(1,2,2);
imagesc(x_v, y_v, mip_filtered'); axis image; colormap hot; colorbar;
title('Laplacian-of-Gaussian Filtered');
xlabel('X (m)'); ylabel('Y (m)');