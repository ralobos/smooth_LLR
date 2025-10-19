% Script to reproduce one of the experiments in  the paper:


clear all;
close all;
clc;

%% Loading data

load('./data/smooth_LLR.mat'); 
% idata_gt - ground truth data
% idata_gt_sc - ground truth single-coil data obtained using a SENSE coil-combination
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

[A, Ah, AhA] = utils.forward_operator(maps, kmask);

%% Creation of patch-extraction operator and its adjoint

patch_size = 8; % patch size for the LLR regularization

[P, Ph] = utils.patch_operators([N1, N2, Nt], patch_size);

%% Gradient function using the hyperbola potential function
% The function also provides the all the terms in the cost function

beta = 5*1e-6; % Regularization parameter

delta = 0.001; % Smoothing parameter for the hyperbola potential function

f_pot = @(x) sqrt(x.^2 + delta^2); % Hyperbola function
f_dot_pot = @(x) x./sqrt(x.^2 + delta^2); % Derivative of the hyperbola function


f_cost_grad = @(x) utils.grad_LLR(x, N1, N2, Nt, ...
                             patch_size, beta, kdata, ...
                             P, Ph, A, Ah, AhA, ...
                             f_pot, f_dot_pot);

%% Huber quadratic majorizer parameters

f_weight_pot = @(x) 1./sqrt(x.^2 + delta^2); % Weight function for the majorizer

step_basis = 0; % Initial step size for the majorizer
niter_maj = 1;  % Number of majorizer iterations
f_maj = @(X, D, step_basis, niter_maj) utils.Huber_quadratic_majorizer(X, D, ...
                                                step_basis, niter_maj, ...
                                                N1, N2, Nt, ...
                                                patch_size, ...
                                                P, ...
                                                A, ...
                                                kdata, ...
                                                beta, ...
                                                f_pot, f_dot_pot, f_weight_pot);

%% Initial reconstruction using data sharing

disp('Computing data-sharing reconstruction...');

cal_length = 8; % Lines in the center of k-space used for calibration
total_lines  = 18; % Total PE lines sampled in each time frame
Raf = 12;

lines_off_acs = total_lines - cal_length;

center = ceil(N2/2)+utils.even_RL(N2);
cal_index = center + [-floor(cal_length/2):floor(cal_length/2)-utils.even_RL(cal_length/2)];
cal_index = cal_index(1:cal_length);

first_idx_acs = cal_index(1);
last_idx_acs = cal_index(end);

kdata_sharing = zeros(N1, N2, Nc, Nt);

kdata_under = kdata;

for t = (Raf + 1 - utils.even_RL(Raf))/2:Nt-(Raf + 1 - utils.even_RL(Raf))/2 
    kdata_sharing(:, :, :, t) = sum(kdata_under(:, :, :, t + [-(Raf-1 + utils.even_RL(Raf))/2 + utils.even_RL(Raf) : (Raf - 1 + utils.even_RL(Raf))/2]), 4); 
    kdata_sharing(:, cal_index, :, t) = kdata_under(:, cal_index, :, t);
end

for t = 1:(Raf + 1 - utils.even_RL(Raf))/2 - 1

    kdata_sharing(:, :, :, t) = sum(kdata_under(:, :, :, t + [-1*(t-1) : Raf - t]), 4); 
    kdata_sharing(:, cal_index, :, t) = kdata_under(:, cal_index, :, t);

end

for t = Nt-(Raf + 1 - utils.even_RL(Raf))/2 + 1 : Nt

    kdata_sharing(:, :, :, t) = sum(kdata_under(:, :, :, t + [-1*(t-1) : Raf - t]), 4); 
    kdata_sharing(:, cal_index, :, t) = kdata_under(:, cal_index, :, t);

end

data_sharing_recon = utils.ift2(kdata_sharing);

% Coil-combination of the data-sharing reconstruction using sensitivity maps

sense_recon_ds = zeros([N1, N2, Nt]);

for t = 1:Nt
    Qlr = sum(conj(maps).*data_sharing_recon(:,:,:,t), 3)./sum(abs(maps).^2, 3);
    Qlr(isnan(Qlr))=0;  
    sense_recon_ds(:,:,t) = Qlr;
end

% NRMSE data-sharing reconstruction

NRMSE_ds = norm(idata_gt_sc(:) - sense_recon_ds(:))/norm(idata_gt_sc(:));

disp(['NRMSE data-sharing reconstruction: ' num2str(NRMSE_ds)]);

% Visualization of data-sharing reconstruction

figure;
imagesc(utils.mdisp(abs(sense_recon_ds))); 
colormap gray;
axis tight;
axis image;
axis off;
title('Data-sharing reconstruction (sum-of-squares) - all frames');

%% Nonlinear conjugate gradient reconstruction

disp(['Starting NCG reconstruction...']);

niter_ncg = 50; % Number of CG iterations

X_pre = sense_recon_ds(:); % Initialization with data-sharing reconstruction

[grad_X_pre, ~, ~, ~, ~] = f_cost_grad(X_pre); % A first gradient calculation

D = -grad_X_pre; % First descending direction

% Vectors to store cost function and NRMSE values

cost_fn_iter = zeros(niter_ncg, 1);
NRMSE_iter = zeros(niter_ncg, 1);

% Figure to visualize reconstruction in each iteration

f_rec_ncg = figure;

for m = 1:niter_ncg

    % ======= Step-size calculation ==============          

    t_ss = tic;
    
    [c0, c1, c2, ~, ~, ~, ~, ~, ~, step_out] = f_maj(X_pre, D, step_basis, niter_maj);

    step_size = step_out - c1/(2*c2);

    disp(['Time step size: ' num2str(toc(t_ss))]);

  % ======= Descending step ==============

    X_next = X_pre + step_size*D;

  % ======= Next descending direction ==============

   [grad_X_next, cost_X_next, ~, ~, ~] = f_cost_grad(X_next);

   beta_next = (norm(grad_X_next)^2)/(norm(grad_X_pre)^2); % Fletcher-Reeves

   D = -grad_X_next + beta_next*D;

   % ======= Update ==============

   error_pre = norm(X_next(:) - X_pre(:))/norm(X_pre(:));
    
   X_pre = X_next;
   grad_X_pre = grad_X_next;

   NRMSE = norm(idata_gt_sc(:) - X_pre(:))/norm(idata_gt_sc(:));

   cost_fn_iter(m) = cost_X_next;
   NRMSE_iter(m) = NRMSE;

   norm_grad = norm(grad_X_pre);

   disp(['Iter: ' int2str(m) ' | Norm grad = ' num2str(norm(grad_X_pre)) ' | step = ' num2str(step_size) ' | cost fn = ' num2str(cost_X_next) ...
         ' | NRMSE = ' num2str(NRMSE) ' | Error Pre = ' num2str(error_pre) ]);

    % Stopping criterion

    if norm_grad < 1e-4 
        break;
    end

    figure(f_rec_ncg); imagesc(utils.mdisp(abs(reshape(X_pre, [N1, N2, Nt])))); 
    axis tight; 
    axis image; 
    axis off;
    colormap gray;
    title(['NCG rec | iter: ' int2str(m)]);

end