% =========================================================================
% function  ph = viewOrthoSlices(img,omega,m,varargin)
%
% Last changed: Lars Ruthotto 2016/02/12
%
% Slice visualization of 3D image.
%
% Input:
%  img   - image data, 3D array
%  omega - image domain
%  m     - voxel size
%
% Output:
%  ph    - plot handle
% =========================================================================
function ph = viewOrthoSlices(img,omega,m,varargin)

phx = []; phy = []; phz = [];
img = reshape(img,m);

roi = [1,m(1),1,m(2),1,m(3)];
sx = round(m(1)/2);
sy = round(m(2)/2);
sz = round(m(3)/2);
cmap = 'jet';
edgeColor = 'none';
for k=1:2:length(varargin),     % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

idx = roi(1):roi(2);
idy = roi(3):roi(4);
idz = roi(5):roi(6);

xn = reshape(getNodalGrid(omega,m),[m+1 3]);
X = xn(:,:,:,1); Y = xn(:,:,:,2); Z = xn(:,:,:,3);
% show x-slice
if not(isempty(sx)),
    getGrid  = @(A,s) squeeze(A(s,idy,idz)+0*A(s+1,idy,idz));
    getSlice = @(A,s) squeeze(A(s,idy,idz));
    phx = surf(getGrid(X,sx), getGrid(Y,sx), getGrid(Z,sx),getSlice(img,sx));
    set(phx,'EdgeColor',edgeColor);
    hold on;
end

% show y-slice
if not(isempty(sy)),
    getGrid  = @(A,s) squeeze(A(idx,s,idz)+0*A(idx,s+1,idz));
    getSlice = @(A,s) squeeze(A(idx,s,idz));
    phy = surf(getGrid(X,sy), getGrid(Y,sy), getGrid(Z,sy),getSlice(img,sy));
    set(phy,'EdgeColor',edgeColor);
    hold on;
end

% show z-slice
if not(isempty(sz)),
    getGrid  = @(A,s) squeeze(A(idx,idy,s)+0*A(idx,idy,s+1));
    getSlice = @(A,s) squeeze(A(idx,idy,s));
    phz = surf(getGrid(X,sz), getGrid(Y,sz), getGrid(Z,sz),getSlice(img,sz));
    set(phz,'EdgeColor',edgeColor);
    hold on;
end

colormap(cmap);
ph = [phx;phy;phz];

function X = getNodalGrid(omega,m)

x1 = []; x2 = []; x3 = [];
h   = (omega(2:2:end)-omega(1:2:end))./m; % voxel size for integration
nu = @(i) (omega(2*i-1)       :h(i):omega(2*i)       )'; % nodal
switch length(omega)/2,
  case 1, x1 = nu(1);
  case 2, [x1,x2] = ndgrid(nu(1),nu(2));
  case 3, [x1,x2,x3] = ndgrid(nu(1),nu(2),nu(3));
end;
X = [x1(:);x2(:);x3(:)];
