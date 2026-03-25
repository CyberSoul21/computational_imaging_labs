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

[NxWall, NyWall, Nt] = size(H);

disp('Size of H:');
disp(size(H));

%% Pre-extract SPAD coordinates
spadX = double(spadPos(:,:,1));
spadY = double(spadPos(:,:,2));
spadZ = double(spadPos(:,:,3));

%% Subsample wall points
stepWall = 8;
sx_list = 1:stepWall:NxWall;
sy_list = 1:stepWall:NyWall;

fprintf('Using %d x %d wall samples\n', length(sx_list), length(sy_list));

%% Small 3D reconstruction grid
Nx = 16;
Ny = 16;
Nz = 16;

x = linspace(volPos(1), volPos(1) + volSize, Nx);
y = linspace(volPos(2), volPos(2) + volSize, Ny);
z = linspace(volPos(3), volPos(3) + volSize, Nz);

G = zeros(Nx, Ny, Nz);

%% Laser wall point
xl = double(reshape(laserPos, [3,1]));

%% Fixed laser device -> laser wall distance
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

                    % voxel -> spad wall
                    d3 = norm(xv - xs);

                    % spad wall -> spad device
                    d4 = norm(xs - spadOrigin);

                    % total optical path
                    tv = d1 + d2 + d3 + d4;

                    % convert to time bin
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

%% Normalize
G = G - min(G(:));
if max(G(:)) > 0
    G = G / max(G(:));
end

%% Visualize slice at the expected object depth
zTarget = 0.5;   % from dataset name Z_d=0.5_...
[~, iz0] = min(abs(z - zTarget));

G_xy_slice = G(:,:,iz0);

figure;
imagesc(x, y, G_xy_slice');
axis image;
axis xy;
colormap hot;
colorbar;
xlabel('x');
ylabel('y');
title(sprintf('XY slice at z = %.3f m', z(iz0)));

%% Also show neighboring slices to inspect depth localization
if iz0 > 1 && iz0 < Nz
    figure;
    subplot(1,3,1);
    imagesc(x, y, G(:,:,iz0-1)');
    axis image; axis xy; colormap hot; colorbar;
    title(sprintf('z = %.3f', z(iz0-1)));

    subplot(1,3,2);
    imagesc(x, y, G(:,:,iz0)');
    axis image; axis xy; colormap hot; colorbar;
    title(sprintf('z = %.3f', z(iz0)));

    subplot(1,3,3);
    imagesc(x, y, G(:,:,iz0+1)');
    axis image; axis xy; colormap hot; colorbar;
    title(sprintf('z = %.3f', z(iz0+1)));
end