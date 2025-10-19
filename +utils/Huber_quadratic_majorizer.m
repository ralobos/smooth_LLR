function [c0_t, c1_t, c2_t, c0_dc, c1_dc, c2_dc, c0_reg, c1_reg, c2_reg, step_out_for_maj] = Huber_quadratic_majorizer(Xin, Pin, step_basis, niter_maj, N1, N2, Nt, patch_size, patches_operator, A, kdata, beta, f_pot, f_dot_pot, f_weight_pot)
% Compute quadratic majorizer coefficients for Huber-based optimization
%
% Input:
%   Xin                 - Vectorized current image estimate
%   Pin                 - Vectorized search direction
%   step_basis          - Initial step size
%   niter_maj           - Number of majorizer iterations
%   N1, N2, Nt         - Image dimensions and time frames
%   patch_size          - Size of local patches
%   patches_operator    - Patch extraction operator
%   A                   - Forward operator
%   kdata               - k-space data
%   beta                - Regularization parameter
%   f_pot, f_dot_pot    - Potential function and derivative
%   f_weight_pot        - Weight function of the potential function
%
% Output:
%   c0_t, c1_t, c2_t    - Total quadratic coefficients: q(step) = c0_t + c1_t*(step-step_ini) + c2_t*(step-step_ini)^2
%   c0_dc, c1_dc, c2_dc - Data consistency term coefficients
%   c0_reg, c1_reg, c2_reg - Regularization term coefficients  
%   step_out_for_maj    - Output step size from majorizer
%
% Note: The 1/2 factor is absorbed into c2_t coefficient
%
% Rodrigo A. Lobos, October 2025


X = reshape(Xin, [N1, N2, Nt]);
P = reshape(Pin, [N1, N2, Nt]);

N1_patch = patch_size^2;
N2_patch = Nt;

q = min(N1_patch, N2_patch);
qm = max(N1_patch, N2_patch);

shifts = -(patch_size/2 - 1) : (patch_size/2); 

% Auxiliary functions

f_diag = @(x) diag(x);

f_c = @(x) repmat(permute(x, [2 1]), [qm, 1]);

step_ini = step_basis;

c2_dc = 0.5*(norm(A(Pin))^2);

for m = 1:niter_maj

    X_basis = X + step_ini*P;
    
    c0_reg = 0;
    c1_reg = 0;
    c2_reg = 0;
    
    for sx = 1:numel(shifts)
        for sy = 1:numel(shifts)
    
            % We shift first 
            
            X_shift = circshift(X_basis, [-shifts(sx) -shifts(sy)]);
            P_shift = circshift(P, [-shifts(sx) -shifts(sy)]);
    
            % Calculation of the three coefficients of the regularizer
            % majorizer
    
            patches_X = patches_operator(X_shift); % dims: [patch_side^2, #time frames, # patches in each image]
            patches_P = patches_operator(P_shift);
    
            N3_patch = size(patches_X, 3);
    
            % SVD calculation of the rectangular matrices
    
            [U, S, V] = pagesvd(patches_X,'vector');
    
            S = squeeze(S);
    
            S_cell = mat2cell(S, q, ones(1,N3_patch)); % Array to cell in order to apply a function to each vector in S
    
            pot_cell = cellfun(f_pot, S_cell , 'UniformOutput', false); % We apply the potential function to each vector in S (in cell format)
    
            weight_pot_cell = cellfun(f_weight_pot, S_cell , 'UniformOutput', false);
    
            dot_pot_cell = cellfun(f_dot_pot, S_cell , 'UniformOutput', false);
    
            G_row = cell2mat(reshape(cellfun(f_c, weight_pot_cell, 'UniformOutput', false), [1 1 N3_patch]));
    
            grad_E_basis = cell2mat(reshape(cellfun(f_diag, dot_pot_cell, 'UniformOutput', false), [1 1 N3_patch]));
     
            M_dir = pagemtimes(pagemtimes(U, 'ctranspose', patches_P, 'none'), 'none', V, 'none');
    
            M_dir_2 = M_dir(1:q, 1:q, :);
            
            c0_reg = c0_reg + (sum(utils.vect(cell2mat(pot_cell))));
            c1_reg = c1_reg + sum(utils.vect(real(conj(M_dir_2).*grad_E_basis)));
            c2_reg = c2_reg + (0.5*sum(utils.vect(G_row.*(abs(M_dir).^2))));

        end
    end

    if m == niter_maj
        step_out_for_maj = step_ini;
    end

    c0_dc = 0.5*norm(A(X_basis(:)) - kdata(:))^2;
    c1_dc = real(A(Pin)'*(A(X_basis(:)) - kdata(:)));
  
    c0_t = c0_dc + beta*c0_reg;
    c1_t = c1_dc + beta*c1_reg;
    c2_t = c2_dc + beta*c2_reg;

    step_ini = step_ini - c1_t/(2*c2_t);
end

end