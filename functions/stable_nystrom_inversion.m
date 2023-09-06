function K = stable_nystrom_inversion(A, B, eps)
% compute A*pinv(B) in a stable manner (A and B denote the AX and Y'AX
% respectively)

[Q1, R] = qr(B, 0);
r = diag(R);
t = (find(abs(r)<eps));
t = min([t; size(R,1)+1]);
R1 = R(1:t-1, :);
[Q2, R2] = qr(R1', 0);
K = A*(Q2/R2')*Q1(:, 1:t-1)';
