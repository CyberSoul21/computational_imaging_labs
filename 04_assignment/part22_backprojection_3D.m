clc;
clear;
close all;

%% Load dataset
dataFolder = '/home/javier/Documents/Computational Imaging/CI_Lab_NLOS_datasets';
fileName = 'Z_d=0.5_l=[1x1]_s=[256x256].mat';
fullPath = fullfile(dataFolder, fileName);

S = load(fullPath);
D = S.data;

H = squeeze(D.data);                    % [NxWall, NyWall, Nt]
spadPos = D.spadPositions;              % [NxWall, NyWall, 3]
laserPos = double(squeeze(D.laserPositions));
laserOrigin = double(D.laserOrigin(:));
spadOrigin = double(D.spadOrigin(:));
volPos = double(D.volumePosition(:));
volSize = double(D.volumeSize);
deltaT = double(D.deltaT);

disp('volumePosition = ');
disp(volPos)

disp('volumeSize = ');
disp(volSize)

disp('Size of H = ');
disp(size(H))

%% Extract SPAD coordinates
spadX = double(spadPos(:,:,1));
spadY = double(spadPos(:,:,2));
spadZ = double(spadPos(:,:,3));

[NxWall, NyWall, Nt] = size(H);

%% Subsample wall points for speed
stepWall = 8;
sx_list = 1:stepWall:NxWall;
sy_list = 1:stepWall:NyWall;

fprintf('Using %d x %d wall samples\n', length(sx_list), length(sy_list));

%% 3D reconstruction grid
Nx = 32;
Ny = 32;
Nz = 11;

% Assume volumePosition is the center of the hidden volume
x = linspace(volPos(1) - volSize/2, volPos(1) + volSize/2, Nx);
y = linspace(volPos(2) - volSize/2, volPos(2) + volSize/2, Ny);

% Narrow slab around expected object depth
zMin = 0.45;
zMax = 0.55;
z = linspace(zMin, zMax, Nz);

G = zeros(Nx, Ny, Nz);

%% Laser wall point
xl = double(reshape(laserPos, [3,1]));

% Fixed distance from laser origin to laser wall point
d1 = norm(laserOrigin - xl);

%% Backprojection
tic;

for ix = 1:Nx
    fprintf('x-slice %d / %d\n', ix, Nx);

    for iy = 1:Ny
        for iz = 1:Nz

            xv = [x(ix); y(iy); z(iz)];
            acc = 0;

            % laser wall -> voxel
            d2 = norm(xl - xv);

            for sx = sx_list
                for sy = sy_list

                    xs = [spadX(sx,sy); spadY(sx,sy); spadZ(sx,sy)];

                    % voxel -> SPAD wall point
                    d3 = norm(xv - xs);

                    % SPAD wall point -> SPAD origin
                    d4 = norm(xs - spadOrigin);

                    % total optical path
                    tv = d1 + d2 + d3 + d4;

                    % convert path length to temporal bin
                    bin = round(tv / deltaT) + 1;

                    if bin >= 1 && bin <= Nt
                        acc = acc + H(sx, sy, bin);
                    end
                end
            end

            G(ix, iy, iz) = acc;
        end
    end
end

elapsedTime = toc;
fprintf('Elapsed time: %.2f seconds\n', elapsedTime);

%% Normalize raw reconstruction
G = G - min(G(:));
if max(G(:)) > 0
    G = G / max(G(:));
end

%% Show XY slice near z = 0.5
zTarget = 0.5;
[~, iz0] = min(abs(z - zTarget));

figure;
imagesc(x, y, G(:,:,iz0)');
axis image;
axis xy;
colormap hot;
colorbar;
xlabel('x');
ylabel('y');
title(sprintf('Raw XY slice at z = %.3f m', z(iz0)));

%% Show neighboring slices
if iz0 > 1 && iz0 < Nz
    figure;
    subplot(1,3,1);
    imagesc(x, y, G(:,:,iz0-1)');
    axis image; axis xy; colormap hot; colorbar;
    xlabel('x'); ylabel('y');
    title(sprintf('z = %.3f', z(iz0-1)));

    subplot(1,3,2);
    imagesc(x, y, G(:,:,iz0)');
    axis image; axis xy; colormap hot; colorbar;
    xlabel('x'); ylabel('y');
    title(sprintf('z = %.3f', z(iz0)));

    subplot(1,3,3);
    imagesc(x, y, G(:,:,iz0+1)');
    axis image; axis xy; colormap hot; colorbar;
    xlabel('x'); ylabel('y');
    title(sprintf('z = %.3f', z(iz0+1)));
end

%% Max projections (optional)
G_xy = squeeze(max(G, [], 3));
G_xz = squeeze(max(G, [], 2));
G_yz = squeeze(max(G, [], 1));

figure;
imagesc(x, y, G_xy');
axis image;
axis xy;
colormap hot;
colorbar;
xlabel('x');
ylabel('y');
title('Max projection XY');

figure;
imagesc(x, z, G_xz');
axis image;
axis xy;
colormap hot;
colorbar;
xlabel('x');
ylabel('z');
title('Max projection XZ');

figure;
imagesc(y, z, G_yz');
axis image;
axis xy;
colormap hot;
colorbar;
xlabel('y');
ylabel('z');
title('Max projection YZ');


%% 3D Laplacian filter
lap3d = zeros(3,3,3);
lap3d(2,2,2) = -6;
lap3d(1,2,2) = 1;
lap3d(3,2,2) = 1;
lap3d(2,1,2) = 1;
lap3d(2,3,2) = 1;
lap3d(2,2,1) = 1;
lap3d(2,2,3) = 1;

G_lap = imfilter(G, -lap3d, 'symmetric');
G_lap(G_lap < 0) = 0;

if max(G_lap(:)) > 0
    G_lap = G_lap / max(G_lap(:));
end

figure;
imagesc(x, y, G_lap(:,:,iz0)');
axis image;
axis xy;
colormap hot;
colorbar;
xlabel('x');
ylabel('y');
title(sprintf('Laplacian filtered XY slice at z = %.3f m', z(iz0)));


%% 3D LoG filtering (Gaussian + Laplacian)

sigma = 0.8;
halfSize = 2;

[xg, yg, zg] = ndgrid(-halfSize:halfSize, -halfSize:halfSize, -halfSize:halfSize);
gauss3d = exp(-(xg.^2 + yg.^2 + zg.^2) / (2*sigma^2));
gauss3d = gauss3d / sum(gauss3d(:));

G_smooth = imfilter(G, gauss3d, 'symmetric');

lap3d = zeros(3,3,3);
lap3d(2,2,2) = -6;
lap3d(1,2,2) = 1;
lap3d(3,2,2) = 1;
lap3d(2,1,2) = 1;
lap3d(2,3,2) = 1;
lap3d(2,2,1) = 1;
lap3d(2,2,3) = 1;

G_log = imfilter(G_smooth, -lap3d, 'symmetric');
G_log(G_log < 0) = 0;

if max(G_log(:)) > 0
    G_log = G_log / max(G_log(:));
end

figure;

subplot(1,2,1);
imagesc(x, y, G(:,:,iz0)');
axis image;
axis xy;
colormap hot;
colorbar;
xlabel('x');
ylabel('y');
title(sprintf('Raw slice at z = %.3f m', z(iz0)));

subplot(1,2,2);
imagesc(x, y, G_log(:,:,iz0)');
axis image;
axis xy;
colormap hot;
colorbar;
xlabel('x');
ylabel('y');
title(sprintf('LoG filtered slice at z = %.3f m', z(iz0)));

%3D representation
title('Volumen 3D')

%volshow(G)
volshow(G, 'RenderingStyle', 'VolumeRendering');