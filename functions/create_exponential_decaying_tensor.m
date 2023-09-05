function T = create_exponential_decaying_tensor(N, rho)

% N: vecotor with the sizes of the tensor 

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
U = cell(1, d);
for i=1:d
    [U{i}, ~] = qr(randn(N(i), n), 0);
end
T  = ttm(S, U); % C times U_k in mode-k for all k.
end
