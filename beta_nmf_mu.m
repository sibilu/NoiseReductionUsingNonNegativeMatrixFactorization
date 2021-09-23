function [D, A, cost] = beta_nmf_mu(S, n_iter, D, A, beta)
for iter = 2:n_iter
% Update D
%D = D .* (( S .* S_ap .^( beta -2) )*A') ./( S_ap .^( beta -1) *A') ;
S_ap = D* A;
% Update A
A = A .* (D'*( S .* S_ap .^( beta -2) )) ./( D'* S_ap .^( beta -1) );
S_ap = D* A;
% Norm -2 normalization
scale = sqrt(sum(D.^2 ,1) );
%D = D .* repmat ( scale.^-1 ,F ,1) ;
A = A .* repmat ( scale',1 , 1549);

% Compute cost value
cost(iter) = betaDivergenceCost (S , S_ap , beta);
end

