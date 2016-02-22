function ERR = TikOPT(lambda, V, bhat, s, V_lbd, xtrue,mask, scale)
%
%    ERR = TikOPT(lambda, V, bhat, s, V_lbd, xtrue, mask, scale)
%
%  This function evaluates the error function for Tikhonov
%  regularization.  
%
%  Input:  lambda -  regularization parameter
%               V - singular vectors of subproblem
%            bhat -  vector U'*b, where U = left singular vectors
%               s -  vector containing the singular values
%           V_lbd - basis vectors from GK bidiagonalization
%           xtrue - true image
%     mask, scale - used to get QSM reconstruction
%
%  Output:  ERR = ||xtrue - V*y_k|| where y_k is the Tikhonov solution to
%  the subproblem
%
%  J.Chung and K. Palmer 3/2014
% Chung and Ruthotto, modified to work for QSM 2016

D = abs(s).^2 + lambda^2;
bhat = conj(s) .* bhat(1:length(s));
xhat = bhat ./ D;
y = V * xhat;
x = V_lbd*y;
ERR = norm(xtrue(:) - mask.*(x/scale));