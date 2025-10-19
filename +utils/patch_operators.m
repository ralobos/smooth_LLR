function [P, Ph] = patch_operators(im_dims, patch_size)
% Create forward and adjoint patch extraction operators
%
% Input:
%   im_dims    - Image dimensions [N1, N2, Nt]
%   patch_size - Size of square patches to extract
%
% Output:
%   P  - Forward operator: extracts patches from image
%   Ph - Adjoint operator: reconstructs image from patches
%
% Note: Assumes image dimensions are divisible by patch_size
%
% Rodrigo A. Lobos, October 2025

N1 = im_dims(1);
N2 = im_dims(2);
Nt = im_dims(3);

f = @(x) reshape(x, [patch_size*patch_size , Nt]);
fh = @(x) reshape(x, [patch_size patch_size Nt]);

number_patches = (N1/patch_size) * (N2/patch_size);

P = @(x) cell2mat(reshape(cellfun(f, mat2cell(x, ...
    repmat(patch_size, [1 N1/patch_size]), repmat(patch_size, [1 N2/patch_size]), [Nt]), ...
    'UniformOutput', false), [1 1 number_patches]));

Ph = @(x) cell2mat(reshape(cellfun(fh, reshape(mat2cell(x, ...
    patch_size*patch_size, Nt, repmat(1, [1 number_patches])), [number_patches 1]), ...
    'UniformOutput', false), [N1/patch_size N2/patch_size]));