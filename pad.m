% =========================================================================
% function  Ip = pad(I,mp)
%
% Last changed: Lars Ruthotto 2016/02/12
%
% Pad image I
%
% Input:
%  I   - image data, 3D array
%  mp  - voxel size of padded volume
%
% Output:
%  Ip  - padded image
% =========================================================================
function Ip = pad(I,mp)
Ip = zeros(mp);
m = size(I);
p = floor((mp-m)/2)+1;
Ip(p(1):p(1)+m(1)-1,p(2):p(2)+m(2)-1,p(3):p(3)+m(3)-1) =I;
end