function B = multilinear_hmt(A, R)

% A input tensor in tensor format, R vector with the desired approximation
% ranks.

N = size(A);
d = length(N);
X = cell(1, d);
Q = cell(1, d);
for i = 1:d
    not_i = N;
    not_i(i) = [];
    %X{i} = randn(prod(not_i), R(i));
    X{i} = randn([not_i, R(i)]);
end
for i = 1:d
    not_i = 1:d;
    not_i(i) = [];
    A_iX_i = tensorprod(A.data, X{i}, not_i, 1:d-1);
    %A_iX_i = tenmat(A, i, not_i).data*X{i};
    [Q{i},~] = qr(A_iX_i, 0); 
end
C = A.data;
for i = 1:d
    C = tensorprod(C, Q{i}, i, 1);
    C = permute(C, [1:i-1, d, i:d-1]);
end
%C = ttm(A, Q, 't'); % A times Q_k^T in mode-k for all k.
B = C;
for i = 1:d
    B = tensorprod(B, Q{i}, i, 2);
    B = permute(B, [1:i-1, d, i:d-1]);
end
%B = ttm(C, Q);