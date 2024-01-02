function exponential_growth_of_the_error(p)
% This function shows that when the error of approximation is big, the
% bounds provided in our paper very descriptive.
% In particular confirming that the presence of an exponential growth of the error is
% possible and not just an artifact of the analysis.
%
% to replicate the experiments in the paper use the following parameters
% for the random generator p: 

rng(p) % 1, 28, 134, 195, 293
%n = 64;
n = 32;
if n == 64
    [U, ~] = qr(randn(50,50), 0);
    U = blkdiag(eye(14), U);
    sigma = diag(10.^(-(0:16/(n-1):16))).*diag(randi([0,1],n,1)*2-1);
elseif n == 32
   [U, ~] = qr(randn(25,25), 0);
    U = blkdiag(eye(7), U);
    sigma = diag(10.^(-(0:16/(n-1):16))).*diag(randi([0,1],n,1)*2-1);
end

rho = 0.3;

%---------------------------
% N = ones(1,2)*n;
% n = min(N);
% d = length(N);
% S = zeros(ones(1,d)*n);
% super_diagonal_indeces = cell(1,d);
% % next part construct a superdiagonal tensor with decaying in the singular
% % values
% for j = 1:n
%     for i = 1:d
%         super_diagonal_indeces{i} = j;
%     end
%     S(super_diagonal_indeces{:}) = rho^j;
% 
% end
% 
% S = tensor(S);
% U_cell = cell(1, d);
% for i=1:d
%     U_cell{i}= U;
% end
% T  = ttm(S, U_cell); % C times U_k in mode-k for all k.
% B = multilinear_nystrom_SRHT(T, [15,15], [5,5], p);
% 
% norma_2 = norm(T-B);

%---------------------------
% N = ones(1,3)*n;
% n = min(N);
% d = length(N);
% S = zeros(ones(1,d)*n);
% super_diagonal_indeces = cell(1,d);
% % next part construct a superdiagonal tensor with decaying in the singular
% % values
% for j = 1:n
%     for i = 1:d
%         super_diagonal_indeces{i} = j;
%     end
%     S(super_diagonal_indeces{:}) = rho^j;
% end
% 
% S = tensor(S);
% U_cell = cell(1, d);
% for i=1:d
%     U_cell{i}= U;
% end
% T  = ttm(S, U_cell); % C times U_k in mode-k for all k.
% 
% B = multilinear_nystrom_SRHT(T, [15,15,15], [5,5,5], p);
% 
% norma_3 = norm(T-B);
%---------------------------
N = ones(1,4)*n;
n = min(N);
d = length(N);
S = zeros(ones(1,d)*n);
super_diagonal_indeces = cell(1,d);
% next part construct a superdiagonal tensor with decaying in the singular
% values
for j = 1:n
    for i = 1:d
        super_diagonal_indeces{i} = j;
    end
    S(super_diagonal_indeces{:}) = rho^j;
end

S = tensor(S);
U_cell = cell(1, d);
for i=1:d
    U_cell{i}= U;
end
T  = ttm(S, U_cell); % C times U_k in mode-k for all k.

B = multilinear_nystrom_SRHT(T, [15,15,15,15], [5,5,5,5], p);
norma_4 = norm(T-B);
%---------------------------------
% N = ones(1,5)*n;
% n = min(N);
% d = length(N);
% S = zeros(ones(1,d)*n);
% super_diagonal_indeces = cell(1,d);
% % next part construct a superdiagonal tensor with decaying in the singular
% % values
% for j = 1:n
%     for i = 1:d
%         super_diagonal_indeces{i} = j;
%     end
%     S(super_diagonal_indeces{:}) = rho^j;
% end
% 
% S = tensor(S);
% U_cell = cell(1, d);
% for i=1:d
%     U_cell{i}= U;
% end
% T  = ttm(S, U_cell); % C times U_k in mode-k for all k.
% 
% B = multilinear_nystrom_SRHT(T, [15,15,15,15,15], [5,5,5,5,5], p);
% norma_5 = norm(T-B);
% 
% 
