clc;
clear;
close all;

dataFolder = '/home/javier/Documents/Computational Imaging/CI_Lab_NLOS_datasets';
fileName = 'Z_d=0.5_l=[1x1]_s=[256x256].mat';

fullPath = fullfile(dataFolder, fileName);
S = load(fullPath);

% Extract the inner struct
D = S.data;

% Extract transient data H
H = D.data;

disp('Original size of H:');
disp(size(H));

% Remove singleton dimensions
H = squeeze(H);

disp('Size of H after squeeze:');
disp(size(H));

% Now H should be [Nx, Ny, Nt]
H_xt = squeeze(sum(H, 2));   % collapse y -> x-t
H_yt = squeeze(sum(H, 1));   % collapse x -> y-t

% Normalize for visualization
if max(H_xt(:)) > 0
    H_xt = H_xt / max(H_xt(:));
end

if max(H_yt(:)) > 0
    H_yt = H_yt / max(H_yt(:));
end

figure;
imagesc(H_xt');
axis xy;
colormap hot;
colorbar;
xlabel('x');
ylabel('time bin');
title(['x-t slice: ', fileName], 'Interpreter', 'none');

figure;
imagesc(H_yt');
axis xy;
colormap hot;
colorbar;
xlabel('y');
ylabel('time bin');
title(['y-t slice: ', fileName], 'Interpreter', 'none');