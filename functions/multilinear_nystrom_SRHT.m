function B = multilinear_nystrom_SRHT(A, R, L, p)

% A input tensor in tensor format, R vector with the desired approximation
% ranks, L vector with the sizes of the oversamples.
N = size(A);
d = length(N);
X = cell(1, d);
Y = cell(1, d);
for i = 1:d
    not_i = N;
    not_i(i) = [];
    X{i} = SRHT_sketching_fixed_rng(prod(not_i), R(i), p);
    Y{i} = SRHT_sketching_fixed_rng(N(i), R(i)+L(i), p);
end
C = ttm(A, Y, 't'); % A times Y_k^T in mode-k for all k.
P = cell(1, d);
Q = cell(1, d);
for i = 1:d
    not_i = 1:d;
    not_i(i) = [];
    A_iX_i = tenmat(A, i, not_i).data*X{i};
    [Q, ~] = qr(A_iX_i, 0);
    err_hmt = norm(A-ttm(A, Q*Q', i));
    err_gn = sqrt(1 + norm(pinv(Y{i}'*Q), 'fro')^2*norm(Y{i}'*(eye(size(Q,1))-Q*Q'), 'fro'));
    P{i} = nystrom_inversion(A_iX_i, Y{i}'*A_iX_i);
    %P{i} = stable_nystrom_inversion(A_iX_i, Y{i}'*A_iX_i, 10*eps);
end 
B = A;
norms = zeros(1,d);
rho = 0.3;
RHO = rho.^(1:1:32);
best_norm = sqrt(sum(RHO(16:end).^2));
for i = 1:d
    B = ttm(B, P{i}*Y{i}', i);
    norms(i) = norm(B-A);
end
    sprintf(['The best approximant has error %d.\n' ...
        'Rho is %d.\n' ...
        'Tau is %d.\n' ...
        'The relative erros are: %d %d %d %e\n' ...
        'Our estimated errors are %d  %d %d  %e.'], ...
        best_norm,err_hmt, err_gn, norms(1)/best_norm, norms(2)/best_norm, norms(3)/best_norm, norms(4)/best_norm, ...
        (err_hmt*(1+err_gn))/best_norm, (err_hmt*(1+err_gn)^2)/best_norm, (err_hmt*(1+err_gn)^3)/best_norm, (err_hmt*(1+err_gn)^4)/best_norm)

B = ttm(C, P);