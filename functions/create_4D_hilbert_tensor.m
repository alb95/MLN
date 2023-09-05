function T = create_4D_hilbert_tensor(n)

% N: vecotor with the sizes of the tensor 

N = ones(1, 4)*n;
% next part construct a superdiagonal tensor with decaying in the singular
% values
T = zeros(N);
for i = 1:n
    for j = 1:n
        for k = 1:n
            for s = 1:n
            T(i,j,k, s) = 1/(i+j+k+s-3);
            end
        end
    end
end
T = tensor(T);