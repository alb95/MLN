function B = sequential_multilinear_hmt(A, R)

% A input tensor in tensor format, R vector with the desired approximation
% ranks, L vector with the sizes of the oversamples.

N = size(A);
d = length(N);
X = cell(1, d);
Q = cell(1, d);
for i = 1:d
    RplusL_prec = R(1:i-1);
    not_i_succ = N(i+1:end);
    X{i} = randn([RplusL_prec, not_i_succ, R(i)]);
    %X{i} = randn(prod(RplusL_prec)*prod(not_i_succ), R(i));
end
C = A.data;
for i = 1:d
    not_i = 1:d;
    not_i(i) = [];
    B_iX_i = tensorprod(C, X{i}, not_i, 1:d-1);
    %B_iX_i = tenmat(C, i, not_i).data*X{i};
    [Q{i},~] = qr(B_iX_i, 0); 
    C = tensorprod(C, Q{i}, i, 1);
    C = permute(C, [1:i-1, d, i:d-1]);
    %C = ttm(C, Q{i}, i, 't');
end 
B = C;
for i = 1:d
    B = tensorprod(B, Q{i}, i, 2);
    B = permute(B, [1:i-1, d, i:d-1]);
end
%B = ttm(C, Q);