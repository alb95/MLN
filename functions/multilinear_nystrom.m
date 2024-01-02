function B = multilinear_nystrom(A, R, L)

% A input tensor in tensor format, R vector with the desired approximation
% ranks, L vector with the sizes of the oversamples.

N = size(A);
d = length(N);
X = cell(1, d);
Y = cell(1, d);
for i = 1:d
    not_i = N;
    not_i(i) = [];
    X{i} = randn([not_i, R(i)]);
    Y{i} = randn([N(i), R(i)+L(i)]);
    %X{i} = randn(prod(not_i), R(i));
    %Y{i} = randn(N(i), R(i)+L(i));
end

C = A.data;
for i = 1:d
    C = tensorprod(C, Y{i}, i, 1);
    C = permute(C, [1:i-1, d, i:d-1]);
end
%C = ttm(A, Y, 't'); % A times Y_k^T in mode-k for all k.
P = cell(1, d);
for i = 1:d
    not_i = 1:d;
    not_i(i) = [];
    A_iX_i = tensorprod(A.data, X{i}, not_i, 1:d-1);
    %A_iX_i = tenmat(A, i, not_i).data*X{i};
    P{i} = nystrom_inversion(A_iX_i, Y{i}'*A_iX_i);
end 
B = C;
for i = 1:d
    B = tensorprod(B, P{i}, i, 2);
    B = permute(B, [1:i-1, d, i:d-1]);
end
%B = ttm(C, P);