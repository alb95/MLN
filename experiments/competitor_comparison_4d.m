% Comparison between:
% 1) HOSVD
% 2) RHOSVD
% 3) MLN (oversample l = floor(0.5r))
% N: size of the tensor, R: multilinear rank of the approximant,
% L: size of oversamples
% sigma: decay_rate for exponential case
% power: decay_rate for polynomial case (= k means sigma(i)=1/(i^k)).     

% experiments parameters:
% 1) N = [100,100,100], power = 1.
% 2) N = [100,100,100], power = 2.
% 3) N = [100,100,100], power = 3.
% 4) N = [100,100,100], sigma = 0.5.
rng(19)
N = [100,100,100, 100];
sigma = 0.5;
%power = 1;
%power = 2;
%power = 3;
%T = create_exponential_decaying_tensor(N, sigma);
%T = create_polynomial_decaying_tensor(N, power);
T = create_4D_hilbert_tensor(100);
Tnorm = norm(T);
E_hosvd = zeros(1,30);
E_rhosvd = zeros(1,30);
E_rsthosvd = zeros(1, 30);
E_mln = zeros(1,30);

for i = 1:30
    r = 2*i;
    B_hosvd = multilinear_svd(T, [r,r,r,r]);
    B_rhosvd = multilinear_hmt(T, [r,r,r,r]);
    B_srhosvd = sequential_multilinear_hmt(T, [r,r,r,r]);
    B_mln = multilinear_nystrom(T, [r,r,r,r], [floor(r/2),floor(r/2), floor(r/2),floor(r/2)]);
    E_hosvd(i) = norm(T-B_hosvd)/Tnorm;
    E_rhosvd(i) = norm(T-B_rhosvd)/Tnorm;
    E_rsthosvd(i) = norm(T-B_srhosvd)/Tnorm;
    E_mln(i) = norm(T-B_mln)/Tnorm;
end

ranks = 2:2:60;
semilogy(ranks, E_hosvd, '-')
hold on
plot(ranks, E_rhosvd,'*')
plot(ranks, E_rsthosvd, '-o')
plot(ranks, E_mln, '--')
legend('HOSVD','RHOSVD','RSTHOSVD','MLN ($\ell=r/2)$', 'Interpreter','latex')
