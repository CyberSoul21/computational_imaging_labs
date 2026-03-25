% =========================================================
%  NLOS Assignment 4 — Complete Solution
%  Sections 2.1, 2.2, 3, and 4
%
%  How to use:
%    1. Set datasetPath to your datasets folder
%    2. Run section by section (use Run Section or Ctrl+Enter)
%       or run the whole file
%    3. Toggle RUN_* flags below to enable/disable sections
% =========================================================
clear; clc; close all;

% ---- USER SETTINGS ----------------------------------------
datasetPath = '/home/javier/Documents/Computational Imaging/CI_Lab_NLOS_datasets/';

RUN_SECTION_21   = true;   % Spatio-temporal slice visualization
RUN_SECTION_22   = true;   % Backprojection + filtering (Z and usaf)
RUN_SECTION_3    = false;%true;   % Confocal bunny reconstruction
RUN_SECTION_4    = false;%true;   % Phasor-field (Morlet wavelet) filtering

% Voxel grid resolution — trade-off: speed vs quality
%   8  → very fast,  coarse
%   16 → fast,       acceptable for testing
%   32 → moderate,  good quality  (recommended)
%   64 → slow,       high quality
Nv   = 32;

% Relay-wall downsampling factor
%   8 → very fast,  1 → full resolution
skip = 8;
% -----------------------------------------------------------


%% ==========================================================
%  SECTION 2.1 — Dataset loading & spatio-temporal slices
% ===========================================================
if RUN_SECTION_21
    fprintf('\n======  SECTION 2.1  ======\n');

    % --- Load the non-time-normalized dataset (Z) ----------
    dsZ = load(fullfile(datasetPath, 'Z_d=0.5_l=[1x1]_s=[256x256].mat'));
    HZ  = dsZ.data.data;          % [Nl_x, Nl_y, Ns_x, Ns_y, Nt]

    % Fix: safely extract scalar fields regardless of vector orientation
    [Nl_xZ, Nl_yZ, Ns_xZ, Ns_yZ, NtZ] = size(HZ);

    % Pick the middle laser point for the slice
    il_mid = round(Nl_xZ/2);
    jl_mid = round(Nl_yZ/2);

    % x-t slice: fix laser to middle, fix SPAD y to middle, vary SPAD x and t
    Hxt_Z = squeeze(HZ(il_mid, jl_mid, :, round(Ns_yZ/2), :));  % [Ns_x, Nt]

    % y-t slice: fix laser to middle, fix SPAD x to middle, vary SPAD y and t
    Hyt_Z = squeeze(HZ(il_mid, jl_mid, round(Ns_xZ/2), :, :));  % [Ns_y, Nt]

    figure('Name','2.1 — Z dataset  (non-time-normalized)','NumberTitle','off',...
           'Position',[50 50 1000 420]);
    subplot(1,2,1);
    imagesc(Hxt_Z); colormap(hot); colorbar; axis tight;
    xlabel('Time bin'); ylabel('SPAD x index');
    title('Z dataset — x–t slice (note skewed wavefront)');

    subplot(1,2,2);
    imagesc(Hyt_Z); colormap(hot); colorbar; axis tight;
    xlabel('Time bin'); ylabel('SPAD y index');
    title('Z dataset — y–t slice');

    % --- Load the time-normalized dataset (usaf) -----------
    dsU = load(fullfile(datasetPath, 'usaf_d=0.5_l=[1x1]_s=[256x256].mat'));
    HU  = dsU.data.data;
    [Nl_xU, Nl_yU, Ns_xU, Ns_yU, NtU] = size(HU);
    il_midU = round(Nl_xU/2);
    jl_midU = round(Nl_yU/2);

    Hxt_U = squeeze(HU(il_midU, jl_midU, :, round(Ns_yU/2), :));
    Hyt_U = squeeze(HU(il_midU, jl_midU, round(Ns_xU/2), :, :));

    figure('Name','2.1 — usaf dataset  (time-normalized)','NumberTitle','off',...
           'Position',[50 520 1000 420]);
    subplot(1,2,1);
    imagesc(Hxt_U); colormap(hot); colorbar; axis tight;
    xlabel('Time bin'); ylabel('SPAD x index');
    title('usaf dataset — x–t slice (wavefront not skewed)');

    subplot(1,2,2);
    imagesc(Hyt_U); colormap(hot); colorbar; axis tight;
    xlabel('Time bin'); ylabel('SPAD y index');
    title('usaf dataset — y–t slice');

    fprintf('Section 2.1 done.\n');
end


%% ==========================================================
%  SECTION 2.2 — Backprojection + filtering
%  Tested on Z dataset and usaf dataset
% ===========================================================
if RUN_SECTION_22

    fprintf('\n======  SECTION 2.2  ======\n');

    datasets_22 = {
        'Z_d=0.5_l=[1x1]_s=[256x256].mat',    false;   % not time-normalized
        'usaf_d=0.5_l=[1x1]_s=[256x256].mat', true;    % time-normalized
    };

    for di = 1:size(datasets_22,1)

        fname      = datasets_22{di,1};
        isTimeNorm = datasets_22{di,2};

        fprintf('\n--- Dataset: %s  (isTimeNorm=%d) ---\n', fname, isTimeNorm);

        ds   = load(fullfile(datasetPath, fname));
        meas = ds.data;

        % --- Unpack -------------------------------------------
        H      = meas.data;
        deltaT = meas.deltaT;
        t0     = meas.t0;

        volPos  = meas.volumePosition(:);   % force column [3x1]
        volSize = meas.volumeSize(:);

        spadPos  = meas.spadPositions;
        laserPos = meas.laserPositions;

        fprintf('  deltaT = %.4f m\n', deltaT);
        fprintf('  Volume centre: [%.3f %.3f %.3f]\n', volPos);
        fprintf('  Volume size:   [%.3f %.3f %.3f]\n', volSize);

        [Nl_x, Nl_y, Ns_x, Ns_y, Nt] = size(H);

        % --- Build voxel grid ---------------------------------

        half = volSize(1) / 2;
        vx = linspace(volPos(1) - half, volPos(1) + half, Nv);
        vy = linspace(volPos(2) - half, volPos(2) + half, Nv);
        vz = linspace(0.45, 0.55, Nv);   % narrow depth slab — much faster

        %vx = linspace(volPos(1)-volSize(1)/2, volPos(1)+volSize(1)/2, Nv);
        %vy = linspace(volPos(2)-volSize(2)/2, volPos(2)+volSize(2)/2, Nv);
        %vz = linspace(volPos(3)-volSize(3)/2, volPos(3)+volSize(3)/2, Nv);

        % --- Reshape relay-wall grids -------------------------
        spadPos_g  = reshape(spadPos,  3, Ns_x, Ns_y);
        laserPos_g = reshape(laserPos, 3, Nl_x, Nl_y);

        % --- Downsample relay wall ----------------------------
        spad_idx  = 1:skip:Ns_x;
        laser_idx = 1:skip:Nl_x;

        H_ds        = H(laser_idx, laser_idx, spad_idx, spad_idx, :);
        spadPos_ds  = spadPos_g(:, spad_idx, spad_idx);
        laserPos_ds = laserPos_g(:, laser_idx, laser_idx);

        [Nl_x_ds, Nl_y_ds, Ns_x_ds, Ns_y_ds, ~] = size(H_ds);
        fprintf('  Downsampled: laser %dx%d  SPAD %dx%d\n', ...
                Nl_x_ds, Nl_y_ds, Ns_x_ds, Ns_y_ds);

        % --- Run backprojection -------------------------------
        G = backproject(H_ds, laserPos_ds, spadPos_ds, ...
                        vx, vy, vz, deltaT, t0, Nt, isTimeNorm, meas);

        % --- Visualize unfiltered projections -----------------
        plot_projections(G, sprintf('Unfiltered — %s  skip=%d  %dx%dx%d', ...
            fname, skip, Nv, Nv, Nv));

        % --- Interactive 3D visualization (MATLAB 2020a+) -----
        figure('Name', sprintf('3D volshow — %s', fname), 'NumberTitle', 'off');
        vs_h = volshow(G, RenderingStyle="MaximumIntensityProjection", Colormap=hot);
        vs_h.Parent.BackgroundColor = [0 0 0];
        vs_h.Parent.GradientColor   = [0 0 0];

        % --- Apply filters ------------------------------------
        % Laplacian
        f_lap = fspecial3('lap');
        G_lap = imfilter(G, -f_lap, 'symmetric');
        G_lap(G_lap < 0) = 0;

        % Laplacian-of-Gaussian (try sigma in [0.5 2.0])
        sigma_log = 1.0;
        f_log     = fspecial3('log', [], sigma_log);
        G_log     = imfilter(G, -f_log, 'symmetric');
        G_log(G_log < 0) = 0;

        % Plot comparison
        plot_filter_comparison(G, G_lap, G_log, ...
            sprintf('Filter comparison — %s  skip=%d  %dx%dx%d', ...
            fname, skip, Nv, Nv, Nv));

    end

    fprintf('Section 2.2 done.\n');
end


%% ==========================================================
%  SECTION 3 — Confocal backprojection (bunny)
% ===========================================================
if RUN_SECTION_3

    fprintf('\n======  SECTION 3  ======\n');

    % --- Non-confocal bunny (5-D H, same as section 2.2) ----
    fname_nc   = 'bunny_d=0.5_l=[1x1]_s=[256x256].mat';
    isTimeNorm_nc = true;   % bunny is time-normalized

    fprintf('Loading non-confocal bunny...\n');
    ds_nc   = load(fullfile(datasetPath, fname_nc));
    meas_nc = ds_nc.data;

    H_nc      = meas_nc.data;
    deltaT_nc = meas_nc.deltaT;
    t0_nc     = meas_nc.t0;
    volPos_nc = meas_nc.volumePosition(:);
    volSize_nc= meas_nc.volumeSize(:);

    [Nl_x, Nl_y, Ns_x, Ns_y, Nt_nc] = size(H_nc);

    vx_nc = linspace(volPos_nc(1)-volSize_nc(1)/2, volPos_nc(1)+volSize_nc(1)/2, Nv);
    vy_nc = linspace(volPos_nc(2)-volSize_nc(2)/2, volPos_nc(2)+volSize_nc(2)/2, Nv);
    vz_nc = linspace(volPos_nc(3)-volSize_nc(3)/2, volPos_nc(3)+volSize_nc(3)/2, Nv);

    spadPos_nc_g  = reshape(meas_nc.spadPositions,  3, Ns_x, Ns_y);
    laserPos_nc_g = reshape(meas_nc.laserPositions, 3, Nl_x, Nl_y);

    spad_idx  = 1:skip:Ns_x;
    laser_idx = 1:skip:Nl_x;

    H_nc_ds        = H_nc(laser_idx, laser_idx, spad_idx, spad_idx, :);
    spadPos_nc_ds  = spadPos_nc_g(:, spad_idx, spad_idx);
    laserPos_nc_ds = laserPos_nc_g(:, laser_idx, laser_idx);

    G_nc = backproject(H_nc_ds, laserPos_nc_ds, spadPos_nc_ds, ...
                       vx_nc, vy_nc, vz_nc, deltaT_nc, t0_nc, Nt_nc, ...
                       isTimeNorm_nc, meas_nc);

    plot_projections(G_nc, sprintf('Non-confocal bunny  skip=%d  %dx%dx%d', skip, Nv,Nv,Nv));

    % Apply LoG filter to non-confocal
    sigma_log = 1.0;
    f_log     = fspecial3('log', [], sigma_log);
    G_nc_log  = imfilter(G_nc, -f_log, 'symmetric');
    G_nc_log(G_nc_log < 0) = 0;
    plot_projections(G_nc_log, sprintf('Non-confocal bunny LoG  skip=%d  %dx%dx%d', skip, Nv,Nv,Nv));

    % --- Confocal bunny (3-D H: [Nc_x, Nc_y, Nt]) -----------
    fname_c   = 'bunny_d=0.5_c=[256x256].mat';
    isTimeNorm_c = true;   % also time-normalized

    fprintf('Loading confocal bunny...\n');
    ds_c   = load(fullfile(datasetPath, fname_c));
    meas_c = ds_c.data;

    H_c      = meas_c.data;       % [Nc_x, Nc_y, Nt]
    deltaT_c = meas_c.deltaT;
    t0_c     = meas_c.t0;
    volPos_c = meas_c.volumePosition(:);
    volSize_c= meas_c.volumeSize(:);

    [Nc_x, Nc_y, Nt_c] = size(H_c);

    vx_c = linspace(volPos_c(1)-volSize_c(1)/2, volPos_c(1)+volSize_c(1)/2, Nv);
    vy_c = linspace(volPos_c(2)-volSize_c(2)/2, volPos_c(2)+volSize_c(2)/2, Nv);
    vz_c = linspace(volPos_c(3)-volSize_c(3)/2, volPos_c(3)+volSize_c(3)/2, Nv);

    % For confocal, spadPositions == laserPositions (co-located)
    confPos_g = reshape(meas_c.spadPositions, 3, Nc_x, Nc_y);

    conf_idx  = 1:skip:Nc_x;
    H_c_ds    = H_c(conf_idx, conf_idx, :);
    confPos_ds= confPos_g(:, conf_idx, conf_idx);

    [Nc_x_ds, Nc_y_ds, ~] = size(H_c_ds);
    fprintf('  Confocal downsampled: %dx%d\n', Nc_x_ds, Nc_y_ds);

    G_c = backproject_confocal(H_c_ds, confPos_ds, vx_c, vy_c, vz_c, ...
                               deltaT_c, t0_c, Nt_c);

    plot_projections(G_c, sprintf('Confocal bunny  skip=%d  %dx%dx%d', skip, Nv,Nv,Nv));

    % Apply LoG filter to confocal
    G_c_log  = imfilter(G_c, -f_log, 'symmetric');
    G_c_log(G_c_log < 0) = 0;
    plot_projections(G_c_log, sprintf('Confocal bunny LoG  skip=%d  %dx%dx%d', skip, Nv,Nv,Nv));

    % --- Side-by-side comparison ------------------------------
    figure('Name','Section 3 — Confocal vs Non-confocal bunny','NumberTitle','off',...
           'Position',[100 100 1200 500]);

    subplot(2,3,1); imagesc(max(G_nc,[],3));     colormap(hot); colorbar; axis image;
    title('Non-confocal  x–y (unfiltered)');
    subplot(2,3,2); imagesc(squeeze(max(G_nc,[],2))); colormap(hot); colorbar; axis image;
    title('Non-confocal  x–z');
    subplot(2,3,3); imagesc(squeeze(max(G_nc,[],1))); colormap(hot); colorbar; axis image;
    title('Non-confocal  y–z');

    subplot(2,3,4); imagesc(max(G_c,[],3));      colormap(hot); colorbar; axis image;
    title('Confocal  x–y (unfiltered)');
    subplot(2,3,5); imagesc(squeeze(max(G_c,[],2)));  colormap(hot); colorbar; axis image;
    title('Confocal  x–z');
    subplot(2,3,6); imagesc(squeeze(max(G_c,[],1)));  colormap(hot); colorbar; axis image;
    title('Confocal  y–z');

    sgtitle('Confocal vs Non-confocal bunny');

    fprintf('Section 3 done.\n');
end


%% ==========================================================
%  SECTION 4 — Phasor-field (Morlet wavelet) filtering
% ===========================================================
if RUN_SECTION_4

    fprintf('\n======  SECTION 4  ======\n');

    % Use the Z dataset (non-time-normalized) as example
    fname_p    = 'Z_d=0.5_l=[1x1]_s=[256x256].mat';
    isTimeNorm_p = false;

    ds_p   = load(fullfile(datasetPath, fname_p));
    meas_p = ds_p.data;

    H_p      = meas_p.data;
    deltaT_p = meas_p.deltaT;
    t0_p     = meas_p.t0;
    volPos_p = meas_p.volumePosition(:);
    volSize_p= meas_p.volumeSize(:);

    [Nl_xp, Nl_yp, Ns_xp, Ns_yp, Nt_p] = size(H_p);

    % Downsample
    spad_idx_p  = 1:skip:Ns_xp;
    laser_idx_p = 1:skip:Nl_xp;

    H_p_ds        = H_p(laser_idx_p, laser_idx_p, spad_idx_p, spad_idx_p, :);
    spadPos_p_g   = reshape(meas_p.spadPositions,  3, Ns_xp, Ns_yp);
    laserPos_p_g  = reshape(meas_p.laserPositions, 3, Nl_xp, Nl_yp);
    spadPos_p_ds  = spadPos_p_g(:, spad_idx_p, spad_idx_p);
    laserPos_p_ds = laserPos_p_g(:, laser_idx_p, laser_idx_p);

    % --- Compute Morlet wavelet parameters -------------------
    % Relay wall SPAD grid spacing (in meters)
    % spadPositions rows: [x;y;z], each column is one point
    % Spacing = distance between adjacent points on the (downsampled) grid
    spadPos_full_g = reshape(meas_p.spadPositions, 3, Ns_xp, Ns_yp);
    dx_spad = abs(spadPos_full_g(1, 2, 1) - spadPos_full_g(1, 1, 1));   % x-spacing
    dy_spad = abs(spadPos_full_g(2, 1, 2) - spadPos_full_g(2, 1, 1));   % y-spacing
    grid_spacing = max(dx_spad, dy_spad) * skip;   % account for downsampling

    lambda_c = 2 * grid_spacing;   % wavelength >= 2 * spacing
    Omega_c  = 1 / lambda_c;       % central frequency
    sigma    = lambda_c;           % sigma in [lambda_c/(2*log2), 2*lambda_c]

    fprintf('  grid_spacing = %.4f m\n', grid_spacing);
    fprintf('  lambda_c     = %.4f m\n', lambda_c);
    fprintf('  Omega_c      = %.4f 1/m\n', Omega_c);
    fprintf('  sigma        = %.4f m\n', sigma);

    % --- Build t array in meter units (optical distance) -----
    t_arr = (0:Nt_p-1) * deltaT_p + t0_p;   % [1 x Nt]

    % --- Morlet wavelet Km(t) --------------------------------
    %   Km(t) = exp(2j*pi*Omega_c*t) * exp(-t^2/(2*sigma^2))
    Km = exp(2j*pi*Omega_c*t_arr) .* exp(-t_arr.^2 / (2*sigma^2));

    % --- Filter H via FFT convolution along temporal dim -----
    % H_p_ds shape: [Nl_x_ds, Nl_y_ds, Ns_x_ds, Ns_y_ds, Nt]
    % Temporal dimension is 5
    temporal_dim = 5;
    Nfft = size(H_p_ds, temporal_dim);

    H_fft  = fft(H_p_ds, Nfft, temporal_dim);
    Km_fft = fft(Km,     Nfft, 2);   % Km is [1 x Nt], fft along dim 2

    % Reshape Km_fft to broadcast across H dimensions [1,1,1,1,Nt]
    Km_fft_bc = reshape(Km_fft, [1, 1, 1, 1, Nfft]);

    H_prime = ifft(H_fft .* Km_fft_bc, Nfft, temporal_dim);

    fprintf('  Morlet filtering done.\n');

    % --- Build voxel grid ------------------------------------
    vx_p = linspace(volPos_p(1)-volSize_p(1)/2, volPos_p(1)+volSize_p(1)/2, Nv);
    vy_p = linspace(volPos_p(2)-volSize_p(2)/2, volPos_p(2)+volSize_p(2)/2, Nv);
    vz_p = linspace(volPos_p(3)-volSize_p(3)/2, volPos_p(3)+volSize_p(3)/2, Nv);

    % --- Run backprojection with H_prime ---------------------
    % H_prime is complex — backproject_phasor accumulates complex values
    G_phasor_complex = backproject_phasor(H_prime, laserPos_p_ds, spadPos_p_ds, ...
                       vx_p, vy_p, vz_p, deltaT_p, t0_p, Nt_p, ...
                       isTimeNorm_p, meas_p);

    % Take absolute value (as specified in the assignment)
    G_phasor = abs(G_phasor_complex);

    plot_projections(G_phasor, sprintf('Phasor (Morlet) — Z  skip=%d  %dx%dx%d', skip, Nv,Nv,Nv));

    % --- Compare: unfiltered vs LoG vs phasor ----------------
    % Re-use unfiltered G from section 2.2 if available, else recompute
    G_plain = backproject(H_p_ds, laserPos_p_ds, spadPos_p_ds, ...
                          vx_p, vy_p, vz_p, deltaT_p, t0_p, Nt_p, ...
                          isTimeNorm_p, meas_p);

    f_log_p  = fspecial3('log', [], 1.0);
    G_log_p  = imfilter(G_plain, -f_log_p, 'symmetric');
    G_log_p(G_log_p < 0) = 0;

    figure('Name','Section 4 — Phasor vs other filters','NumberTitle','off',...
           'Position',[100 100 1400 400]);
    vols_cmp = {G_plain, G_log_p, G_phasor};
    ttls_cmp = {'Unfiltered', 'LoG filter', 'Phasor (Morlet)'};
    for c = 1:3
        subplot(1,3,c);
        imagesc(max(vols_cmp{c},[],3)); colormap(hot); colorbar; axis image;
        title([ttls_cmp{c} '  x–y']);
    end
    sgtitle(sprintf('Filter comparison (Z dataset, skip=%d, %dx%dx%d)', skip, Nv,Nv,Nv));

    fprintf('Section 4 done.\n');
end

fprintf('\n===== ALL DONE =====\n');


%% ==========================================================
%  LOCAL FUNCTIONS
% ===========================================================

% ----------------------------------------------------------
%  backproject  —  non-confocal backprojection (5-D H)
%    H_ds         [Nl_x_ds, Nl_y_ds, Ns_x_ds, Ns_y_ds, Nt]
%    laserPos_ds  [3, Nl_x_ds, Nl_y_ds]
%    spadPos_ds   [3, Ns_x_ds, Ns_y_ds]
%    vx,vy,vz     1-D coordinate arrays  (length Nv each)
%    deltaT       temporal bin width (meters)
%    t0           time offset of first bin (meters)
%    Nt           number of temporal bins
%    isTimeNorm   bool: true = d1 and d4 already subtracted
%    meas         full dataset struct (needed for laserOrigin/spadOrigin)
% ----------------------------------------------------------
function G = backproject(H_ds, laserPos_ds, spadPos_ds, ...
                         vx, vy, vz, deltaT, t0, Nt, isTimeNorm, meas)

    Nv     = numel(vx);
    G      = zeros(Nv, Nv, Nv);
    [Nl_x_ds, Nl_y_ds, Ns_x_ds, Ns_y_ds, ~] = size(H_ds);

    tic;
    fprintf('  Running backprojection...\n');

    for ix = 1:Nv
        for iy = 1:Nv
            for iz = 1:Nv

                xv  = [vx(ix); vy(iy); vz(iz)];
                acc = 0;

                for il_x = 1:Nl_x_ds
                    for il_y = 1:Nl_y_ds

                        xl = laserPos_ds(:, il_x, il_y);
                        d2 = norm(xv - xl);

                        if ~isTimeNorm
                            d1 = norm(xl - meas.laserOrigin(:));
                        else
                            d1 = 0;
                        end

                        for is_x = 1:Ns_x_ds
                            for is_y = 1:Ns_y_ds

                                xs = spadPos_ds(:, is_x, is_y);
                                d3 = norm(xv - xs);

                                if ~isTimeNorm
                                    d4 = norm(xs - meas.spadOrigin(:));
                                else
                                    d4 = 0;
                                end

                                tv  = d1 + d2 + d3 + d4;
                                bin = round((tv - t0) / deltaT) + 1;

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
        fprintf('    slice ix = %d / %d\n', ix, Nv);
    end

    elapsed = toc;
    fprintf('  Backprojection done in %.1f s.\n', elapsed);
end


% ----------------------------------------------------------
%  backproject_confocal  —  confocal backprojection (3-D H)
%    H_ds         [Nc_x_ds, Nc_y_ds, Nt]
%    confPos_ds   [3, Nc_x_ds, Nc_y_ds]  (xl == xs for every point)
%    vx,vy,vz     1-D coordinate arrays
%    deltaT, t0   temporal parameters (meters)
%    Nt           number of temporal bins
%  Note: confocal datasets are time-normalized → no d1/d4
% ----------------------------------------------------------
function G = backproject_confocal(H_ds, confPos_ds, vx, vy, vz, deltaT, t0, Nt)

    Nv = numel(vx);
    G  = zeros(Nv, Nv, Nv);
    [Nc_x_ds, Nc_y_ds, ~] = size(H_ds);

    tic;
    fprintf('  Running CONFOCAL backprojection...\n');

    for ix = 1:Nv
        for iy = 1:Nv
            for iz = 1:Nv

                xv  = [vx(ix); vy(iy); vz(iz)];
                acc = 0;

                for ic_x = 1:Nc_x_ds
                    for ic_y = 1:Nc_y_ds

                        xc = confPos_ds(:, ic_x, ic_y);

                        % Confocal: xl == xs, so tv = d2 + d3
                        d2 = norm(xv - xc);   % laser wall → voxel
                        d3 = norm(xv - xc);   % voxel → SPAD wall (same point)
                        tv = d2 + d3;         % = 2 * norm(xv - xc)

                        bin = round((tv - t0) / deltaT) + 1;

                        if bin >= 1 && bin <= Nt
                            acc = acc + H_ds(ic_x, ic_y, bin);
                        end

                    end
                end

                G(ix, iy, iz) = acc;
            end
        end
        fprintf('    slice ix = %d / %d\n', ix, Nv);
    end

    elapsed = toc;
    fprintf('  Confocal backprojection done in %.1f s.\n', elapsed);
end


% ----------------------------------------------------------
%  backproject_phasor  —  same as backproject but accepts
%  complex H_ds (result of Morlet filtering)
%  Returns complex G — caller should take abs()
% ----------------------------------------------------------
function G = backproject_phasor(H_ds, laserPos_ds, spadPos_ds, ...
                                vx, vy, vz, deltaT, t0, Nt, isTimeNorm, meas)

    Nv     = numel(vx);
    G      = complex(zeros(Nv, Nv, Nv));
    [Nl_x_ds, Nl_y_ds, Ns_x_ds, Ns_y_ds, ~] = size(H_ds);

    tic;
    fprintf('  Running PHASOR backprojection (complex H)...\n');

    for ix = 1:Nv
        for iy = 1:Nv
            for iz = 1:Nv

                xv  = [vx(ix); vy(iy); vz(iz)];
                acc = complex(0);

                for il_x = 1:Nl_x_ds
                    for il_y = 1:Nl_y_ds

                        xl = laserPos_ds(:, il_x, il_y);
                        d2 = norm(xv - xl);

                        if ~isTimeNorm
                            d1 = norm(xl - meas.laserOrigin(:));
                        else
                            d1 = 0;
                        end

                        for is_x = 1:Ns_x_ds
                            for is_y = 1:Ns_y_ds

                                xs = spadPos_ds(:, is_x, is_y);
                                d3 = norm(xv - xs);

                                if ~isTimeNorm
                                    d4 = norm(xs - meas.spadOrigin(:));
                                else
                                    d4 = 0;
                                end

                                tv  = d1 + d2 + d3 + d4;
                                bin = round((tv - t0) / deltaT) + 1;

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
        fprintf('    slice ix = %d / %d\n', ix, Nv);
    end

    elapsed = toc;
    fprintf('  Phasor backprojection done in %.1f s.\n', elapsed);
end


% ----------------------------------------------------------
%  plot_projections  —  3 max-intensity projections
% ----------------------------------------------------------
function plot_projections(G, fig_title)
    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 1200 380]);

    subplot(1,3,1);
    imagesc(max(G, [], 3)); colormap(hot); colorbar; axis image;
    title('x–y  (front)'); xlabel('y'); ylabel('x');

    subplot(1,3,2);
    imagesc(squeeze(max(G, [], 2))); colormap(hot); colorbar; axis image;
    title('x–z  (side)'); xlabel('z'); ylabel('x');

    subplot(1,3,3);
    imagesc(squeeze(max(G, [], 1))); colormap(hot); colorbar; axis image;
    title('y–z  (top)'); xlabel('z'); ylabel('y');

    sgtitle(fig_title, 'Interpreter', 'none');
end


% ----------------------------------------------------------
%  plot_filter_comparison  —  3x3 grid (filters x views)
% ----------------------------------------------------------
function plot_filter_comparison(G, G_lap, G_log, fig_title)
    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 1200 750]);

    titles_rows = {'Unfiltered', 'Laplacian', 'LoG'};
    vols        = {G, G_lap, G_log};

    for r = 1:3
        V = vols{r};

        subplot(3,3,(r-1)*3+1);
        imagesc(max(V,[],3)); colormap(hot); colorbar; axis image;
        title([titles_rows{r} '  x–y']); xlabel('y'); ylabel('x');

        subplot(3,3,(r-1)*3+2);
        imagesc(squeeze(max(V,[],2))); colormap(hot); colorbar; axis image;
        title([titles_rows{r} '  x–z']); xlabel('z'); ylabel('x');

        subplot(3,3,(r-1)*3+3);
        imagesc(squeeze(max(V,[],1))); colormap(hot); colorbar; axis image;
        title([titles_rows{r} '  y–z']); xlabel('z'); ylabel('y');
    end

    sgtitle(fig_title, 'Interpreter', 'none');
end