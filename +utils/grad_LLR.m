function [grad_fn, cost_fn, cost_dc, cost_reg, cost_lambda_reg] = grad_LLR(Xin, N1, N2, Nt, patch_size, lambda, kdata, patches_operator, patches_adj_operator, A_operator, Ah_operator, AhA_operator, f_pot, f_dot_pot)

X = reshape(Xin, [N1, N2, Nt]);

cost_reg = 0;
grad_reg = zeros(prod(size(X)),1);

N1_patch = patch_size^2;
N2_patch = Nt;

q = min(N1_patch, N2_patch);

f_diag = @(x) diag(x);

shifts = -(patch_size/2 - 1) : (patch_size/2); 

for sx = 1:numel(shifts)
    for sy = 1:numel(shifts)
        
        X_shift = circshift(X, [-shifts(sx) -shifts(sy)]);

        patches_X = patches_operator(X_shift);

        N3_patch = size(patches_X, 3);
        
        [U, S, V] = pagesvd(patches_X, 'econ', 'vector');

        S = squeeze(S);

        S_cell = mat2cell(S, q, ones(1,N3_patch));

        pot_cell = cellfun(f_pot, S_cell , 'UniformOutput', false);

        dot_pot_cell = cellfun(f_dot_pot, S_cell , 'UniformOutput', false);

        pot_dot_S = cell2mat(reshape(cellfun(f_diag, dot_pot_cell, 'UniformOutput', false), [1 1 N3_patch]));
        
        cost_reg = cost_reg + sum(utils.vect(cell2mat(pot_cell)));
        
        UpotS = pagemtimes(U,pot_dot_S);
        grad_patches = pagemtimes(UpotS,"none",V,"ctranspose");

        grad_reg = grad_reg + utils.vect(circshift(patches_adj_operator(grad_patches), [shifts(sx) shifts(sy)]));
    end
end

grad_dc = AhA_operator(utils.vect(X)) - Ah_operator(kdata(:));

grad_fn = grad_dc + lambda*grad_reg;

cost_dc = 0.5*norm(A_operator(utils.vect(X)) - kdata(:))^2;

cost_lambda_reg = lambda*cost_reg;

cost_fn = cost_dc + cost_lambda_reg;

end