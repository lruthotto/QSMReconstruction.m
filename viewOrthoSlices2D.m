% =========================================================================
% function  ph = viewOrthoSlices2D(img,omega,m,varargin)
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
function viewOrthoSlices2D(img,omega,m,varargin)

if nargin==0
    runMinimalExample;
    return;
end

slices = ceil(m./2);
color = [.99 0.99 .99];

h = (omega(2:2:end)-omega(1:2:end))./m;
for k=1:2:length(varargin),     % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;
% create big image
img = reshape(img,m);
Ixz = squeeze(img(:,slices(2),:)); % xz view
Iyz = squeeze(img(slices(1),:,:)); % yz view
Ixy = squeeze(img(:,:,slices(3))); % xy view
I   = [Ixz' flipud(Iyz)'; flipud(Ixy') nan(m(2),m(2))];

% put together x and y axis of plot
x = omega(1):h(1):omega(2);
y = omega(3):h(2):omega(4);
z = omega(5):h(3):omega(6);
xa = [x(:);y(2:end)'+omega(2)];
ya = [z(:);y(2:end)'+omega(6);];

Ip = zeros(size(I,1)+1,size(I,2)+1);
Ip(1:end-1,1:end-1)=I;
[Xa,Ya] = meshgrid(xa,ya);
ph = pcolor(Xa,Ya,Ip);
set(ph, 'EdgeColor', 'none');
hold on;
plot([0;max(xa)],[omega(6);omega(6)],'LineWidth',5,'Color',color);
plot([omega(2); omega(2)],[0;max(ya)],'LineWidth',5,'Color',color);
plot([0;max(xa)],[z(slices(3));z(slices(3))],'--','LineWidth',2,'Color',color);
plot([0;max(x)],[y(slices(2));y(slices(2))]+omega(6),'--','LineWidth',2,'Color',color);
plot([x(slices(1));x(slices(1))],[0;max(ya)],'--','LineWidth',2,'Color',color);
plot([y(slices(2));y(slices(2))]+omega(2),[0;max(z)],'--','LineWidth',2,'Color',color);


function runMinimalExample

m      = [8 14 24];
omega  = [0 2 0 4 0 1];
h      = (omega(2:2:end)-omega(1:2:end))./m;

x = h(1)/2:h(1):omega(2);
y = h(2)/2:h(2):omega(4);
z = h(3)/2:h(3):omega(6);

[X,Y,Z] = ndgrid(x,y,z);

figure(1); clf;
viewOrthoSlices2D(flipdim(X,1),omega,m);


