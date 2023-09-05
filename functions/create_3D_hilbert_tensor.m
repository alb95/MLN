function T = create_3D_hilbert_tensor(n)

% N: vecotor with the sizes of the tensor 

N = ones(1, 3)*n;
% next part construct a superdiagonal tensor with decaying in the singular
% values
T = zeros(N);
for i = 1:n
    for j = 1:n
        for k = 1:n
            T(i,j,k) = 1/(i+j+k-2);
        end
    end
end
T = tensor(T);