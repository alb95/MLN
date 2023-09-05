function K = nystrom_inversion(A, B)
% compute A*pinv(B) in a stable manner (A and B denote the AX and Y'AX
% respectively)

[Q, R] = qr(B, 0);
K = (A/R)*Q';
