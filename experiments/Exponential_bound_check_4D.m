%n = 64;
n = 32;
if n == 64
    [U, ~] = qr(randn(50,50), 0);
    U = blkdiag(eye(14), U);
    sigma = diag(10.^(-(0:16/(n-1):16))).*diag(randi([0,1],n,1)*2-1);
elseif n == 32
   [U, ~] = qr(randn(25,25), 0);
    U = blkdiag(eye(7), U);
    sigma = diag(10.^(-(0:16/(n-1):16))).*diag(randi([0,1],n,1)*2-1);
end

%[U, ~] = qr(randn(32,32),0);

% N: vecotor with the sizes of the tensor 

rho = 0.5;
N = ones(1,4)*n;
n = min(N);
d = length(N);
S = zeros(ones(1,d)*n);
super_diagonal_indeces = cell(1,d);
% next part construct a superdiagonal tensor with decaying in the singular
% values
for j = 1:n
    for i = 1:d
        super_diagonal_indeces{i} = j;
    end
    S(super_diagonal_indeces{:}) = rho^j;
end

S = tensor(S);
U_cell = cell(1, d);
for i=1:d
    U_cell{i}= U;
end
T  = ttm(S, U_cell); % C times U_k in mode-k for all k.

C = multilinear_nystrom(T, [15,15,15,15], [5,5,5,5]);
B = multilinear_nystrom_SRHT(T, [15,15,15,15], [5,5,5,5]);
%B = multilinear_nystrom_DCT(T, [15,15,15,15], [5,5,5,5]);
