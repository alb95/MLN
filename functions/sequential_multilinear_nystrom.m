function B = sequential_multilinear_nystrom(A, R, L)

% A input tensor in tensor format, R vector with the desired approximation
% ranks, L vector with the sizes of the oversamples.

%eps = delta*1e-16*norm(A);
N = size(A);
d = length(N);
X = cell(1, d);
Y = cell(1, d);
for i = 1:d
    RplusL_prec = R(1:i-1)+L(1:i-1);
    not_i_succ = N(i+1:end);
    X{i} = randn([RplusL_prec, not_i_succ, R(i)]);
    Y{i} = randn([N(i), R(i)+L(i)]);
    %X{i} = randn(prod(RplusL_prec)*prod(not_i_succ), R(i));
    %Y{i} = randn(N(i), R(i)+L(i));
end
P = cell(1, d);
C = A.data;
for i = 1:d
    not_i = 1:d;
    not_i(i) = [];
    B_iX_i = tensorprod(C, X{i}, not_i, 1:d-1);
    %B_iX_i = tenmat(C, i, not_i).data*X{i};
    P{i} =  nystrom_inversion(B_iX_i, Y{i}'*B_iX_i);
    C = tensorprod(C, Y{i}, i, 1);
    C = permute(C, [1:i-1, d, i:d-1]);
    %C = ttm(C, Y{i}, i, 't');
end 
B = C;
for i = 1:d
    B = tensorprod(B, P{i}, i, 2);
    B = permute(B, [1:i-1, d, i:d-1]);
end
%B = ttm(C, P);