% Comparison of the errors of approximation of MLN 
% for different dimensions d, but same decay in the singular values:
% 1) d = 2
% 2) d = 3
% 3) d = 4
% 4) d = 5

% N = ones(1,d)*n : size of the tensor, R: multilinear rank of the approximant,
% L: size of oversamples
% sigma: decay_rate for exponential case
% power: decay_rate for polynomial case (= k means sigma(i)=1/(i^k)).     

% experiments parameters:
% 1) n = 35, power = 3.
% 2) n = 35, sigma = 0.5.

n = 35;
power = 3;
sigma = 0.5;

T_2 = create_polynomial_decaying_tensor([n,n], power);
T_3 = create_polynomial_decaying_tensor([n,n,n], power);
T_4 = create_polynomial_decaying_tensor([n,n,n,n], power);
T_5 = create_polynomial_decaying_tensor([n,n,n,n,n], power);

%T_2 = create_exponential_decaying_tensor([n,n], sigma);
%T_3 = create_exponential_decaying_tensor([n,n,n], sigma);
%T_4 = create_exponential_decaying_tensor([n,n,n,n], sigma);
%T_5 = create_exponential_decaying_tensor([n,n,n,n,n], sigma);

E_2 = zeros(1,n-5); 
E_3 = zeros(1,n-5); 
E_4 = zeros(1,n-5); 
E_5 = zeros(1,n-5); 

for i = 1:n-5
    r = i;
    B_2 = multilinear_nystrom(T_2, [r,r], ones(1,2)*3);
    B_3 = multilinear_nystrom(T_3, [r,r,r], ones(1,3)*3);
    B_4 = multilinear_nystrom(T_4, [r,r,r, r], ones(1,4)*3);
    B_5 = multilinear_nystrom(T_5, [r,r,r,r,r], ones(1,5)*3);
    E_2(i) = norm(T_2-B_2);
    E_3(i) = norm(T_3-B_3);
    E_4(i) = norm(T_4-B_4);
    E_5(i) = norm(T_5-B_5);
end

ranks = 1:n-5;
semilogy(ranks, E_2, '-')
hold on
plot(ranks, E_3, '.-')
plot(ranks, E_4, '-o')
plot(ranks, E_5,'-*')
legend('$d=2$','$d=3$', '$d=4$', '$d=5$', 'Interpreter','latex')

