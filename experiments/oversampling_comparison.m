% Comparison between different oversampling sizes:
% 1) l = 0
% 2) l = r+5
% 3) l = 0.5r
% 4) l = r
% N: size of the tensor, R: multilinear rank of the approximant,
% L: size of oversamples
% sigma: decay_rate for exponential case
% power: decay_rate for polynomial case (= k means sigma(i)=1/(i^k)).     
% size(T) = [100,100,100], decay rate = 0.7.

% on a 3d-tensor T with exponential decay in the singular values.
% N: size of the tensor, R: multilinear rank of the approximant,
% sigma: decay_rate, L: size of oversamples
% experiments parameter:
% 1) N = [100,100,100], sigma = 0.7.
% 2) N = [100,100,100], power = 1.

N = [100,100,100];
sigma = 0.7;
power = 1;
T = create_exponential_decaying_tensor(N, sigma);
%T = create_polynomial_decaying_tensor(N, power);

E_hosvd = zeros(1,25); 
E_0 = zeros(1,25);
E_2 = zeros(1,25);
E_half = zeros(1,25);
E_r = zeros(1,25);

for i = 1:25
    r = 4*i;
    R = [r,r,r];
    B_hosvd = multilinear_svd(T, R);
    B_0 = multilinear_nystrom(T, R, [0,0,0]);
    B_2 = multilinear_nystrom(T, R, [2,2,2]);
    B_half = multilinear_nystrom(T,R, [floor(r/2),floor(r/2), floor(r/2)]);
    B_r = multilinear_nystrom(T, R, R);
    E_hosvd(i) = norm(T-B_hosvd);
    E_0(i) = norm(T-B_0);
    E_2(i) = norm(T-B_2);
    E_half(i) = norm(T-B_half);
    E_r(i) = norm(T-B_r);
end

ranks = 4:4:100;
semilogy(ranks, E_0, '-')
hold on
plot(ranks, E_2, 'o')
plot(ranks, E_half,'*')
plot(ranks, E_r, '--')
plot(ranks, E_hosvd, '.-')
legend('$\ell=0$','$\ell=2$', '$\ell=r/2$', '$\ell=r$','HOSVD',  'Interpreter','latex')
