% =========================================================================
% function  I = unpad(Ip,m)
%
% Last changed: Lars Ruthotto 2016/02/12
%
% Remove padding from image Ip
%
% Input:
%  Ip   - image data, 3D array
%  m    - voxel size of unpadded volume
%
% Output:
%  I    - unpadded image
% =========================================================================
function I = unpad(Ip,m)
mp = size(Ip);
p = floor((mp-m)/2)+1;
I = Ip(p(1):p(1)+m(1)-1,p(2):p(2)+m(2)-1,p(3):p(3)+m(3)-1) ;

