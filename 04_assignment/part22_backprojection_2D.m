clc;
clear;
close all;

%% Load dataset
dataFolder = '/home/javier/Documents/Computational Imaging/CI_Lab_NLOS_datasets';
fileName = 'Z_d=0.5_l=[1x1]_s=[256x256].mat';
fullPath = fullfile(dataFolder, fileName);

S = load(fullPath);
D = S.data;

H = squeeze(D.data);                    % [256, 256, 4048]
spadPos = D.spadPositions;              % [256, 256, 3]
laserPos = double(squeeze(D.laserPositions));
laserOrigin = double(D.laserOrigin(:));
spadOrigin = double(D.spadOrigin(:));
volPos = double(D.volumePosition(:));
volSize = double(D.volumeSize);
deltaT = double(D.deltaT);

disp('Size of H:');
disp(size(H));

%% Pre-extract SPAD coordinates
spadX = double(spadPos(:,:,1));
spadY = double(spadPos(:,:,2));
spadZ = double(spadPos(:,:,3));

%% Dimensions
[NxWall, NyWall, Nt] = size(H);

%% Subsample wall points for speed
stepWall = 8;
sx_list = 1:stepWall:NxWall;
sy_list = 1:stepWall:NyWall;

fprintf('Using %d x %d wall samples\n', length(sx_list), length(sy_list));

%% Small reconstruction grid
Nx = 16;
Ny = 16;

x = linspace(volPos(1), volPos(1) + volSize, Nx);
y = linspace(volPos(2), volPos(2) + volSize, Ny);
z0 = volPos(3) + 0.5 * volSize;

G = zeros(Nx, Ny);

%% Laser wall point
xl = double(reshape(laserPos, [3,1]));

%% Fixed distance from laser device to laser wall
d1 = norm(laserOrigin - xl);

%% Backprojection
tic;

for ix = 1:Nx
    fprintf('Row %d / %d\n', ix, Nx);

    for iy = 1:Ny
        xv = [x(ix); y(iy); z0];
        acc = 0;

        % laser wall point -> voxel
        d2 = norm(xl - xv);

        for sx = sx_list
            for sy = sy_list
                xs = [spadX(sx,sy); spadY(sx,sy); spadZ(sx,sy)];

                % voxel -> spad wall point
                d3 = norm(xv - xs);

                % spad wall point -> spad device
                d4 = norm(xs - spadOrigin);

                % total path length
                tv = d1 + d2 + d3 + d4;

                % convert to time bin
                bin = round(tv / deltaT) + 1;

                if bin >= 1 && bin <= Nt
                    acc = acc + H(sx, sy, bin);
                end
            end
        end

        G(ix, iy) = acc;
    end
end

elapsedTime = toc;
fprintf('Elapsed time: %.2f seconds\n', elapsedTime);

%% Normalize
if max(G(:)) > 0
    G = G / max(G(:));
end

%% Show reconstruction
figure;
imagesc(x, y, G');
axis image;
axis xy;
colormap hot;
colorbar;
xlabel('x');
ylabel('y');
title('2D backprojection reconstruction');