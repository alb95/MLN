% Accuracy comparison between MLN and SMLN:
% 1) MLN
% 2) SMLN
% on a 3d-tensor T with exponential decay in the singular values.
% N: size of the tensor, R: multilinear rank of the approximant,
% sigma: decay_rate, L: size of oversamples
% experiments parameter:
% 1) N = [70,70,70], sigma = 0.1, L = [0,0,0].
% 2) N = [70,70,70], sigma = 0.1, L = [3,3,3].

N = [70,70,70];
sigma = 0.1;
L = [0,0,0];
T = create_exponential_decaying_tensor(N, sigma);
E_HOSVD = zeros(1,N(1));
E_MLN = zeros(1,N(1)); 
E_SMLN_1 = zeros(1,N(1));
E_SMLN_10 = zeros(1,N(1));
for i = 1:70
    r = i;
    R = [r,r,r];
    B_HOSVD = multilinear_svd(T, R);
    B_MLN = multilinear_nystrom(T, R, L);
    B_SMLN_1 = stabilized_multilinear_nystrom(T, R, L, 1);
    B_SMLN_10 = stabilized_multilinear_nystrom(T, R, L, 10);
    E_HOSVD(i) = norm(T-B_HOSVD);
    E_MLN(i) = norm(T-B_MLN);
    E_SMLN_1(i) = norm(T-B_SMLN_1);
    E_SMLN_10(i) = norm(T-B_SMLN_10);
end

ranks = 1:N(1);
semilogy(ranks, E_HOSVD, '-')
hold on
plot(ranks, E_MLN, '-o')
plot(ranks, E_SMLN_1, '.-')
plot(ranks, E_SMLN_10, '-*')
legend('HOSVD', 'MLN', 'SMLN-1', 'SMLN-10')
