function B = multilinear_svd(A, R)
% HOSVD with ranks given by R

% A input tensor in tensor format, R vector with the desired approximation
% ranks.
N = size(A);
d = length(N);
Q = cell(1, d);
for i = 1:d
    A_i = tenmat(A, i).data;
    [Q{i},~,~] = svd(A_i, 'econ');
    Q{i} = Q{i}(:, 1:R(i));
end
C = ttm(A, Q, 't'); % A times Q_k^T in mode-k for all k.
B = ttm(C, Q);