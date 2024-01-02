% This script show that MLN with structured sketchings (Khatri-Rao form), 
% recompresses faster the sum of two tensors in Tucker form. 

% N: vecotor with the sizes of the tensor 
% rho1, rho2: exponential decay rates
% R: vector representing the multilinear rank of the tensor
% L: oversampling parameter
% ATTENTION: for the moment, only for tensor with same sizes and ranks


N = [3000,3000,3000];
R = [60,60,60];
L = [10,10,10];
rho1 = 0.4;
rho2 = 0.5;

% CONSTRUCT THE TENSOR TO BE SUMMED.
d = length(N);
n = min(N);
r = min(R);
l = min(L);
S1 = zeros(ones(1,d)*r);
S2 = zeros(ones(1,d)*r);
super_diagonal_indeces = cell(1,d);
% next part construct a superdiagonal tensor with decaying in the singular
% values
for j = 1:r
    for i = 1:d
        super_diagonal_indeces{i} = j;
    end
    S1(super_diagonal_indeces{:}) = rho1^j;
    S2(super_diagonal_indeces{:}) = rho2^j;
end

U1 = cell(1, d);
U2 = cell(1, d);
for i=1:d
    [U1{i}, ~] = qr(randn(N(i), r), 0);
    [U2{i}, ~] = qr(randn(N(i), r), 0);
end

S = zeros(ones(1,d)*2*r);
index1 = cell(1,d);
index2 = cell(1,d);
for i = 1:d
    index1{i} = 1:r;
end
for i = 1:d
    index2{i} = r+1:2*r;
end
S(index1{:}) = S1;
S(index2{:}) = S2;

U = cell(1, d);
for i=1:d
    U{i} = [U1{i}, U2{i}];
end

%---------------------------------
% RHOSVD
tic

X = cell(1, d);
Wx = cell(1, d); % U{i}'*X{i}

for i = 1:d 
    X{i} = randn(n, r);
    Wx{i} = U{i}'*X{i};
end

Psi = cell(1, d);
for i = 1:d 
    not_i = 1:d;
    not_i(i) = [];
    Psi{i} = khatrirao(Wx{not_i});
    Psi{i} = reshape(Psi{i}, [ones(1,d-1)*(2*r), r]);
end

Q = cell(1, d);
for i = 1:d
    not_i = 1:d;
    not_i(i) = [];
    S_iPsi_i = tensorprod(S, Psi{i}, not_i, 1:d-1);
    U_iS_iPsi_i = U{i}*S_iPsi_i; 
    [Q{i}, ~] = qr(U_iS_iPsi_i, 0);
end 

for i = 1:d 
    K{i} = Q{i}'*U{i};
end

D = S;
for i = 1:d
    D = tensorprod(D, K{i}', i, 1);
    D = permute(D, [1:i-1, d, i:d-1]);
end

t_rhosvd = toc;

%---------------------------------

% RSTHOSVD TODO
% tic
% 
% X = cell(1, d);
% Wx = cell(1, d); % U{i}'*X{i}
% 
% for i = 1:d 
%     X{i} = randn(n, r);
%     Wx{i} = U{i}'*X{i};
% end
% 
% Psi = cell(1, d);
% for i = 1:d 
%     not_i = 1:d;
%     not_i(i) = [];
%     Psi{i} = khatrirao(Wx{not_i});
%     Psi{i} = reshape(Psi{i}, [ones(1,d-1)*(2*r), r]);
% end
% 
% Q = cell(1, d);
% for i = 1:d
%     not_i = 1:d;
%     not_i(i) = [];
%     S_iPsi_i = tensorprod(S, Psi{i}, not_i, 1:d-1);
%     U_iS_iPsi_i = U{i}*S_iPsi_i; 
%     [Q{i}, ~] = qr(U_iS_iPsi_i, 0);
% end 
% 
% for i = 1:d 
%     K{i} = Q{i}'*U{i};
% end
% 
% D = S;
% for i = 1:d
%     D = tensorprod(D, K{i}', i, 1);
%     D = permute(D, [1:i-1, d, i:d-1]);
% end
% 
% t_rsthosvd = toc;

%---------------------------------

% MLN

tic

X = cell(1, d);
Wx = cell(1, d); % U{i}'*X{i}

for i = 1:d 
    X{i} = randn(n, r+l);
    Wx{i} = U{i}'*X{i};
end

Psi = cell(1, d);
for i = 1:d 
    not_i = 1:d;
    not_i(i) = [];
    Psi{i} = khatrirao(Wx{not_i});
    Psi{i} = reshape(Psi{i}, [ones(1,d-1)*(2*r), r+l]);
end

Y = cell(1,d);
Wy = cell(1, d); % U{i}'*Y{i}
for i = 1:d 
    Y{i} = randn(n, r);
    Wy{i} = U{i}'*Y{i};
end

C = S;
for i = 1:d
    C = tensorprod(C, Wy{i}, i, 1);
    C = permute(C, [1:i-1, d, i:d-1]);
end

P = cell(1, d);
for i = 1:d
    not_i = 1:d;
    not_i(i) = [];
    S_iPsi_i = tensorprod(S, Psi{i}, not_i, 1:d-1);
    P{i} = nystrom_inversion(S_iPsi_i, Wy{i}'*S_iPsi_i);
    P{i} = U{i}*P{i};
end 

t_mln = toc;
%------------------------------

% CONSTRUCTION OF THE FULL TENSOR

% classic sum
%for i = 1:d
%    S = tensorprod(S, U{i}, i, 2);
%    S = permute(S, [1:i-1, d, i:d-1]);
%end

% RHOSVD sum
%for i = 1:d
%    D = tensorprod(D, Q{i}, i, 2);
%    D = permute(D, [1:i-1, d, i:d-1]);
%end

% MLN sum
%B = C;
%for i = 1:d
%    B = tensorprod(B, P{i}, i, 2);
%    B = permute(B, [1:i-1, d, i:d-1]);
%end
