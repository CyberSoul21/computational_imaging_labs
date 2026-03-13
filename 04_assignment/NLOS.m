clear; clc; close all;
D = load("CI_Lab_NLOS_datasets\Z_d=0.5_l=[1x1]_s=[256x256].mat");

fn = fieldnames(D);
dataset = D.(fn{1});

% extraer histograma
H = squeeze(dataset.data);

[nx,ny,nt] = size(H)

% slices centrales
x_mid = round(nx/2);
y_mid = round(ny/2);

Hx = squeeze(H(:,y_mid,:));
Hy = squeeze(H(x_mid,:,:));

% mejorar visibilidad
Hx = log(Hx + 1);
Hy = log(Hy + 1);

figure
imagesc(Hx)
axis xy
colormap hot
colorbar
title('x-t slice')

figure
imagesc(Hy)
axis xy
colormap hot
colorbar
title('y-t slice')