function W = SRHT_sketching(n, k)
H0 = [1, 1; 1, -1]/sqrt(2);
H = H0;
for i = 1:log2(n)-1
    H = kron(H, H0);
end
D = diag((rand(1,n)<.5)*2 - 1);
W = sqrt(n/k)*D*H(:, randperm(n,k));