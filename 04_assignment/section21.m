%% =========================================================
%  NLOS Assignment 4 — Section 2.1  (corrected)
% ==========================================================
clear; clc; close all;

datasetPath = '/home/javier/Documents/Computational Imaging/CI_Lab_NLOS_datasets/';

files = {
    'Z_d=0.5_l=[1x1]_s=[256x256].mat',
    'usaf_d=0.5_l=[1x1]_s=[256x256].mat',
};
labels = {'Z (not time-normalized)', 'USAF (time-normalized)'};

for i = 1:length(files)

    fprintf('\n--- Loading: %s ---\n', files{i});
    ds = load(fullfile(datasetPath, files{i}));

    % The struct is nested: ds.data contains all fields
    meas = ds.data;

    % H shape: [Nl_x, Nl_y, Ns_x, Ns_y, Nt]
    %   Nl_x, Nl_y = laser grid  (1x1 here, single laser point)
    %   Ns_x, Ns_y = SPAD grid   (256x256 detectors)
    %   Nt          = temporal bins
    H = meas.data;
    fprintf('Raw H size: %s\n', mat2str(size(H)));

    % Squeeze out the singleton laser dimensions → [Ns_x, Ns_y, Nt]
    H_sq = squeeze(H);   % [256, 256, 4048]

    Ns_x = size(H_sq, 1);
    Ns_y = size(H_sq, 2);
    Nt   = size(H_sq, 3);
    fprintf('Squeezed H: Ns_x=%d  Ns_y=%d  Nt=%d\n', Ns_x, Ns_y, Nt);

    % --- Spatio-temporal slices at the centre of each spatial axis ---
    mid_x = round(Ns_x / 2);
    mid_y = round(Ns_y / 2);

    % x-t slice: for every x, take the central y row  → [Ns_x, Nt]
    slice_xt = squeeze(H_sq(:, mid_y, :));

    % y-t slice: for every y, take the central x column → [Ns_y, Nt]
    slice_yt = squeeze(H_sq(mid_x, :, :));

    % --- Plot ---
    figure('Name', labels{i}, 'NumberTitle', 'off', ...
           'Position', [100, 100, 1200, 450]);

    subplot(1,2,1);
    imagesc(slice_xt);
    colormap(hot); colorbar;
    xlabel('Temporal bin'); ylabel('SPAD x-index');
    title([labels{i} '  |  x–t slice']);
    axis tight;

    subplot(1,2,2);
    imagesc(slice_yt);
    colormap(hot); colorbar;
    xlabel('Temporal bin'); ylabel('SPAD y-index');
    title([labels{i} '  |  y–t slice']);
    axis tight;

    sgtitle(labels{i}, 'FontSize', 13, 'FontWeight', 'bold');
end

fprintf('\nDone.\n');