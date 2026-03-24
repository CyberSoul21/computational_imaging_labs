%% =========================================================
%  NLOS Assignment 4 — Section 2.2
%  Backprojection reconstruction + filtering
% ==========================================================
clear; clc; close all;

datasetPath = '/home/javier/Documents/Computational Imaging/CI_Lab_NLOS_datasets/';

%% ---------------------------------------------------------
%  Choose which dataset to reconstruct.
%  Z  → NOT time-normalized  (include d1 + d4)
%  usaf → IS  time-normalized  (skip d1 + d4)
% ----------------------------------------------------------
datasetFile   = 'Z_d=0.5_l=[1x1]_s=[256x256].mat';
isTimeNorm    = false;   % <-- set true for usaf dataset

% Load
ds   = load(fullfile(datasetPath, datasetFile));
meas = ds.data;

%% ---------------------------------------------------------
%  1. Unpack measurement struct
% ----------------------------------------------------------
H      = meas.data;          % [Nl_x, Nl_y, Ns_x, Ns_y, Nt]
deltaT = meas.deltaT;        % temporal bin width in METERS (optical path)
t0     = meas.t0;            % time offset of first bin (usually 0)

% Relay-wall grid positions (3 x N matrices, each column is [x;y;z])
spadPos  = meas.spadPositions;   % 3 x (Ns_x*Ns_y)
laserPos = meas.laserPositions;  % 3 x (Nl_x*Nl_y)

% Hidden volume extent
volPos  = meas.volumePosition;   % [x;y;z] centre of hidden volume
volSize = meas.volumeSize;       % [sx;sy;sz] full extent

fprintf('deltaT = %.4f m\n', deltaT);
fprintf('Volume centre: [%.2f, %.2f, %.2f] m\n', volPos);
fprintf('Volume size:   [%.2f, %.2f, %.2f] m\n', volSize);

%% ---------------------------------------------------------
%  2. Build the voxel grid for the hidden scene
% ----------------------------------------------------------
% Resolution — start small to check correctness, then increase
Nv = 32;   % voxels per side  (try 16, 32, 64)

% Voxel coordinates span [volPos - volSize/2 , volPos + volSize/2]
vx = linspace(volPos(1) - volSize(1)/2,  volPos(1) + volSize(1)/2,  Nv);
vy = linspace(volPos(2) - volSize(2)/2,  volPos(2) + volSize(2)/2,  Nv);
vz = linspace(volPos(3) - volSize(3)/2,  volPos(3) + volSize(3)/2,  Nv);

fprintf('Voxel grid: %d x %d x %d\n', Nv, Nv, Nv);

%% ---------------------------------------------------------
%  3. Unpack H and prepare relay-wall grids
% ----------------------------------------------------------
% H shape: [Nl_x, Nl_y, Ns_x, Ns_y, Nt]
[Nl_x, Nl_y, Ns_x, Ns_y, Nt] = size(H);

% Reshape positions to grids for easy indexing
% spadPos / laserPos are 3 x N (each column is one point [x;y;z])
% Reshape to 3 x Ns_x x Ns_y and 3 x Nl_x x Nl_y
spadPos_g  = reshape(spadPos,  3, Ns_x, Ns_y);   % [3, Ns_x, Ns_y]
laserPos_g = reshape(laserPos, 3, Nl_x, Nl_y);   % [3, Nl_x, Nl_y]

%% ---------------------------------------------------------
%  4. (Optional) Downsample relay wall to speed up testing
% ----------------------------------------------------------
% Use every `skip`-th point. Set skip=1 for full resolution.
skip = 8;   % try 8 for quick test, 4 for better quality, 1 for full

spad_idx  = 1:skip:Ns_x;
laser_idx = 1:skip:Nl_x;

H_ds        = H(laser_idx, laser_idx, spad_idx, spad_idx, :);
spadPos_ds  = spadPos_g(:, spad_idx, spad_idx);
laserPos_ds = laserPos_g(:, laser_idx, laser_idx);

[Nl_x_ds, Nl_y_ds, Ns_x_ds, Ns_y_ds, ~] = size(H_ds);
fprintf('Downsampled grid: laser %dx%d, SPAD %dx%d\n', ...
        Nl_x_ds, Nl_y_ds, Ns_x_ds, Ns_y_ds);

%% ---------------------------------------------------------
%  5. Backprojection
%
%  For every voxel xv:
%    G(xv) = sum over all (xl, xs) of  H(xl, xs, tv)
%
%  where tv = d1 + d2 + d3 + d4   (full path in meters)
%    or tv = d2 + d3               (time-normalized datasets)
%
%  The temporal bin index is:  bin = round((tv - t0) / deltaT) + 1
% ----------------------------------------------------------
G = zeros(Nv, Nv, Nv);   % reconstruction volume

tic;
fprintf('Running backprojection...\n');

for ix = 1:Nv
for iy = 1:Nv
for iz = 1:Nv

    xv = [vx(ix); vy(iy); vz(iz)];   % current voxel position [3x1]
    acc = 0;                           % accumulator for this voxel

    for il_x = 1:Nl_x_ds
    for il_y = 1:Nl_y_ds

        xl = laserPos_ds(:, il_x, il_y);   % laser point [3x1]

        % d2: distance from laser wall point to voxel
        d2 = norm(xv - xl);

        % d1: distance from laser device to xl
        % For 1x1 laser grids, laserOrigin gives the device position
        if ~isTimeNorm
            d1 = norm(xl - meas.laserOrigin);
        else
            d1 = 0;
        end

        for is_x = 1:Ns_x_ds
        for is_y = 1:Ns_y_ds

            xs = spadPos_ds(:, is_x, is_y);   % SPAD point [3x1]

            % d3: distance from voxel to SPAD wall point
            d3 = norm(xv - xs);

            % d4: distance from xs to SPAD device
            if ~isTimeNorm
                d4 = norm(xs - meas.spadOrigin);
            else
                d4 = 0;
            end

            % Total optical path length (in meters)
            tv = d1 + d2 + d3 + d4;

            % Convert to bin index (deltaT is already in meters)
            bin = round((tv - t0) / deltaT) + 1;

            % Accumulate if bin is within the valid range
            if bin >= 1 && bin <= Nt
                acc = acc + H_ds(il_x, il_y, is_x, is_y, bin);
            end

        end
        end
    end
    end

    G(ix, iy, iz) = acc;

end
end
    % Progress report every row
    fprintf('  slice ix = %d / %d\n', ix, Nv);
end

elapsed = toc;
fprintf('Backprojection done in %.1f seconds.\n', elapsed);

%% ---------------------------------------------------------
%  6. Visualize — 2D max-intensity projections
% ----------------------------------------------------------
% Maximum intensity projection along depth (z) → front view (x-y)
G_xy = max(G, [], 3);   % max over z
G_xz = max(G, [], 2);   % max over y  (side view)
G_yz = max(G, [], 1);   % max over x  (top view)
G_yz = squeeze(G_yz);

figure('Name', 'Backprojection — unfiltered', 'NumberTitle', 'off', ...
       'Position', [100, 100, 1200, 380]);

subplot(1,3,1);
imagesc(G_xy); colormap(hot); colorbar; axis image;
title('x–y projection (front)'); xlabel('y'); ylabel('x');

subplot(1,3,2);
imagesc(squeeze(G_xz)); colormap(hot); colorbar; axis image;
title('x–z projection (side)'); xlabel('z'); ylabel('x');

subplot(1,3,3);
imagesc(G_yz); colormap(hot); colorbar; axis image;
title('y–z projection (top)'); xlabel('z'); ylabel('y');

sgtitle(sprintf('Unfiltered backprojection  |  %dx%dx%d voxels  |  skip=%d  |  %.0fs', ...
        Nv, Nv, Nv, skip, elapsed));

%% ---------------------------------------------------------
%  7. Filter — Laplacian and Laplacian-of-Gaussian
% ----------------------------------------------------------

% --- Laplacian filter ---
f_lap   = fspecial3('lap');
G_lap   = imfilter(G, -f_lap, 'symmetric');
G_lap(G_lap < 0) = 0;

% --- Laplacian-of-Gaussian filter ---
% sigma controls the smoothing; try values between 0.5 and 2.0
sigma_log = 1.0;
f_log     = fspecial3('log', [], sigma_log);
G_log     = imfilter(G, -f_log, 'symmetric');
G_log(G_log < 0) = 0;

% --- Plot filtered projections ---
figure('Name', 'Filtered reconstructions', 'NumberTitle', 'off', ...
       'Position', [100, 100, 1200, 750]);

titles_rows = {'Unfiltered', 'Laplacian', 'Laplacian-of-Gaussian'};
vols        = {G, G_lap, G_log};

for r = 1:3
    V = vols{r};
    subplot(3, 3, (r-1)*3 + 1);
    imagesc(max(V, [], 3)); colormap(hot); colorbar; axis image;
    title([titles_rows{r} '  x–y']); xlabel('y'); ylabel('x');

    subplot(3, 3, (r-1)*3 + 2);
    imagesc(squeeze(max(V, [], 2))); colormap(hot); colorbar; axis image;
    title([titles_rows{r} '  x–z']); xlabel('z'); ylabel('x');

    subplot(3, 3, (r-1)*3 + 3);
    imagesc(squeeze(max(V, [], 1))); colormap(hot); colorbar; axis image;
    title([titles_rows{r} '  y–z']); xlabel('z'); ylabel('y');
end

sgtitle(sprintf('Filter comparison  |  %dx%dx%d voxels  |  skip=%d', ...
        Nv, Nv, Nv, skip));

fprintf('All done.\n');