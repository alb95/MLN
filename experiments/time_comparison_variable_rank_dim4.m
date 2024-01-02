% Time Comparison between:
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
% 4) N = [100,100,100], sigma = 0.7.

N = [70,70, 70, 70];
sigma = 0.7;
power = 1;
% power = 2;
% power = 3;
%T = create_exponential_decaying_tensor(N, sigma);
T = create_polynomial_decaying_tensor(N, power);

E_hosvd = zeros(1,20);
E_rhosvd = zeros(1,20);
E_mln = zeros(1,20);
t_hosvd = zeros(1,20);
t_rhosvd = zeros(1,20);
t_seq_rhosvd = zeros(1,20);
t_seq_mln = zeros(1,20);
t_mln = zeros(1,20);
    for i = 1:20
        r = 2*i;
        %tic
        %B_hosvd = multilinear_svd(T, [r,r,r,r]);
        %t_hosvd(i) = toc;
        t_rhosvd(i) = timeit(@() multilinear_hmt(T, [r,r,r,r]));
        t_mln(i) = timeit(@() multilinear_nystrom(T, [r,r,r,r], [floor(r/2),floor(r/2), floor(r/2), floor(r/2)]));
        t_seq_rhosvd(i) = timeit(@() sequential_multilinear_hmt(T, [r,r,r,r]));
        %t_seq_mln(i) = timeit(@() sequential_multilinear_nystrom(T, [r,r,r,r], [floor(r/2),floor(r/2), floor(r/2),floor(r/2)]));
        %E_hosvd(i) = norm(T-B_hosvd);
        %E_rhosvd(i) = norm(T-B_rhosvd);
        %E_mln(i) = norm(T-B_mln);
    end

%ranks = 3:3:60;
%semilogy(ranks, E_rhosvd, '-')
%hold on
%plot(ranks, E_mln, '--')
%legend('HOSVD','RHOSVD', 'MLN ($\ell=r/2)$', 'Interpreter','latex')
%legend('RHOSVD', 'MLN ($\ell=r/2)$', 'Interpreter','latex')
%title('a')

figure
ranks = 2:2:40;
semilogy(ranks, t_rhosvd, '-')
hold on
plot(ranks, t_seq_rhosvd, '--')
plot(ranks, t_mln, '-*')
%legend('HOSVD','RHOSVD', 'MLN ($\ell=r/2)$', 'Interpreter','latex')
legend('RHOSVD', 'SRHOSVD', 'MLN ($\ell=r/2)$', 'Interpreter','latex')
title('time comparison')


