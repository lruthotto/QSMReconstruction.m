% =========================================================================
% function  V = sdiag(v)
%
% Last changed: Lars Ruthotto 2016/02/12
%
% shortcut for generating sparse diagonal matrix from vector
%
% Input:
%  v   - vector
%
% Output:
%  V  - sparse diagonal matrix
% =========================================================================
function V = sdiag(v)
    V = spdiags(v(:),0,numel(v),numel(v));
end