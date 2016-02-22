function DP = TikDP(lambda, V, bhat, s, nLevel, B, beta, N)
%
%    DP = TikDP(lambda, V, bhat, s, nLevel, B, beta, N)
%
%  This function evaluates the DP function for Tikhonov
%  regularization.
%
%  Input:  lambda -  regularization parameter
%            bhat -  vector U'*b, where U = left singular vectors
%               s -  vector containing the singular values
%          nLevel -  standard deviation of the noise
%                   If nLevel = 'est'=, then we estimate the noise level
%               B - bidiagonal matrix
%            beta - defines right hand side vector
%               N - length of noise vector
%                   For problems where the norm of the residual is included
%                   as nLevel, we ignore the input N
%
%  Output:     DP - abs(||b-Ax_lam||^2- tau*nLevel^2 * N)
%
%      Note: nLevel^2 * N is the expected value of ||res||_2^2
%             tau is a safety parameter
%
%  J.Chung and K. Palmer 3/2014

tau = 1;

bhat = conj(s) .* bhat(1:length(s));
vector = beta*eye(size(B,1),1);
for i = 1:length(lambda)
  
  D = abs(s).^2 + lambda(i)^2;
  xhat = bhat ./ D;
  x = V * xhat;
  
  res = vector - B*x;
  DP(i) = abs(norm(res)^2- tau * (nLevel^2));
end

