function W = SRHT_sketching_fixed_rng(n, k, p)
rng(p)
H0 = [1; 1]/sqrt(2);
H1 = [1; -1]/sqrt(2);
HR = zeros(n,k);
for j = 1:k
    x = randi([0,1]);
    if x == 0
        H = H0;
    else 
        H = H1;
    end
    for i = 1:log2(n)-1
        x = randi([0,1]);
        if x == 0
            H =  kron(H, H0);
        else
            H =  kron(H, H1);
        end
    end
    HR(:, j) = H; 
end
D = spdiags((rand(n,1)<.5)*2 - 1, 0, n, n);

W = sqrt(n/k)*D*HR;