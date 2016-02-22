%% Fig21_PSFs.m
%  Compares the point spread function for image deblurring and QSM 
%   reconstruction example.  Orthogonal slices of the convolution kernels
%   are provided in both the spatial and frequency domains.
%  
%   This script generates Figure 2.1 in
%   "Computational Methods for Image Reconstruction" - Chung and Ruthotto
%       NMR Biomedicine Special Issue: QSM, 2016 
%
% Chung & Ruthotto (2016)

fontSize = 10;

matrix_size = [256   256   146];
voxel_size = [0.9375 0.9375 1];
B0_dir = [0.0665919 0 0.997780000000000];

omega = zeros(1,6);
omega(1:2:end) = -matrix_size./2;
omega(2:2:end) = matrix_size./2;

%% PSF for deblurring
[PSF, center] = Gauss( [5 5 5] , matrix_size);

fig = figure(21); clf;
fig.Name = 'Fig 2.1, Chung and Ruthotto (2016): PSF and Deconvolution Kernel';
subplot(2,3,1)
viewOrthoSlices(PSF(129-3:129+4,129-3:129+4,74-3:74+4),.5*[-8 8 -8 8 -8 8],.5*[16 16 16]);
xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'FontSize',fontSize)
xlim = caxis();
axis tight
title('deblurring, PSF')

view([146,32])
colorbar('XTickLabel',{sprintf('%1.e',min(xlim(1))); sprintf('%1.e',(max(xlim(2))))},'XTick',[(min(xlim(1))); (max(xlim(2)))])
% In frequency domain
S = fftn(circshift(PSF,1-center)); % Singular values
S_shift = fftshift(S); % Shift so that (0,0,0) is the center of the freq domain

subplot(2,3,2);
viewOrthoSlices(abs(S_shift),omega,matrix_size);
xlabel('k_x')
ylabel('k_y')
zlabel('k_z')
colorbar
set(gca,'FontSize',fontSize)
axis tight
xlim = caxis();
colorbar('XTickLabel',{sprintf('%1.1f',min(xlim(1))); sprintf('%1.1f',(max(xlim(2))))},'XTick',[(min(xlim(1))); (max(xlim(2)))])
title('deblurring, PSF in frequency domain')

subplot(2,3,3)
viewOrthoSlices(log10(abs(S_shift)),omega,matrix_size);
xlabel('k_x')
ylabel('k_y')
zlabel('k_z')
colorbar
set(gca,'FontSize',fontSize)
axis tight
xlim = caxis();
colorbar('XTickLabel',{sprintf('%1.1f',min(xlim(1))); sprintf('%1.1f',(max(xlim(2))))},'XTick',[(min(xlim(1))); (max(xlim(2)))])
title('deblurring, log of Fourier coefficients')


%% PSF for dipole kernel (1)
% (1) generate the e-values directly
D = dipole_kernel(matrix_size, voxel_size, B0_dir);
D_shift = fftshift(D);

% Invert to get the PSF in image space
PSF_D = ifftn(D);
PSF_D_shift = fftshift(PSF_D);

subplot(2,3,4)
viewOrthoSlices(real(PSF_D_shift(129-3:129+4,129-3:129+4,74-3:74+4)),[-4 4 -4 4 -4 4],[8 8 8]);
xlabel('x')
ylabel('y')
zlabel('z')
colorbar
view([146,32])
set(gca,'FontSize',fontSize)
axis tight
xlim = caxis();
colorbar('XTickLabel',{sprintf('%1.1f',min(xlim(1))); sprintf('%1.1f',(max(xlim(2))))},'XTick',[(min(xlim(1))); (max(xlim(2)))])
title('QSM, PSF')

subplot(2,3,5)
viewOrthoSlices(D_shift,omega,matrix_size);
xlabel('k_x')
ylabel('k_y')
zlabel('k_z')
colorbar
set(gca,'FontSize',fontSize)
axis tight
xlim = caxis();
colorbar('XTickLabel',{sprintf('%1.1f',min(xlim(1))); sprintf('%1.1f',(max(xlim(2))))},'XTick',[(min(xlim(1))); (max(xlim(2)))])
title('QSM, PSF in frequency domain')

subplot(2,3,6)
viewOrthoSlices(log10(abs(D_shift)),omega,matrix_size);
xlabel('k_x')
ylabel('k_y')
zlabel('k_z')
set(gca,'FontSize',fontSize)
axis tight
caxis([-3.5 0]);
xlim = caxis();
colorbar('XTickLabel',{sprintf('<%1.1f',min(xlim)); sprintf('%1.1f',(max(xlim)))},'XTick',[(min(xlim)); (max(xlim))])
title('QSM, log of Fourier coefficients')

