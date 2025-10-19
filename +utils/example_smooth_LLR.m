% Script to reproduce one of the experiments in  the paper:


clear all;
close all;
clc;

%% Loading data

load('./data/smooth_LLR.mat'); 
% idata_gt - ground truth data
% kmask    - sampling mask
% maps     - sensitivity maps

%% Data dimensions

[N1, N2, Nc, Nt] = size(idata_gt);  % N1 x N2 : image dimensions
                                    % Nc      : number of coils
                                    % Nt      : number of time frames

%% sum-of-squares coil combination

idata_gt_sos = squeeze(sqrt(sum(abs(idata_gt).^2,3)));

%% Visualization of all the frames

figure;
imagesc(utils.mdisp(abs(idata_gt_sos))); 
colormap gray; 
axis tight;
axis image;
axis off;
title('Ground truth data - all frames');

%% Visualization of k-space mask (all frames)

figure;
imagesc(utils.mdisp(squeeze(kmask(:, :, 1, :)))); 
colormap gray; 
axis tight;
axis image;
axis off;
title('k-space sampling mask - all frames');

%% Ground-truth data in k-space

kdata_gt = utils.ft2(idata_gt);

%% Retrospectively undersampled k-space data

kdata = kdata_gt .* kmask;

%% sum-of-squares coil combination of undersampled data

idata_under_sos = squeeze(sqrt(sum(abs(utils.ift2(kdata)).^2,3)));

figure;
imagesc(utils.mdisp(abs(idata_under_sos))); 
colormap gray;
axis tight;
axis image;
axis off;
title('Undersampled data - all frames');

%% Creation of the forward operator, its adjoint, and the composition of both (A, Ah, AhA)

