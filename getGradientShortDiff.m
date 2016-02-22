% =========================================================================
% function [L,D1,D2,D3] = getGradientShortDiff(m,h)
%
% Last changed: Lars Ruthotto 2016/02/12
%
% Get short difference discretization of 3D gradient operator.
%
% Input:
%  m   - number of voxels
%  h   - voxel size
%
% Output:
%  G   = [D1; D2; D3] image gradient
%  Di  - ith partial derivative operator, i=1,2,3
%
% =========================================================================
function [G,D1,D2,D3] = getGradientShortDiff(m,h)

dx = @(i) spdiags(ones(m(i),1)*[-1 1], [0 1], m(i)-1,m(i))/h(i);

D1 = kron(speye(m(3)*m(2)),dx(1));
D2 = kron(speye(m(3)),kron(dx(2),speye(m(1))));
D3 = kron(dx(3), speye(m(2)*m(1)));
G  = [D1; D2; D3];

