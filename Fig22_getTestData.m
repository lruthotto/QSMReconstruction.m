%% Fig22_getTestData.m
%  Generates an image deblurring example and a QSM reconstruction example.
%  
%   This script generates data in Figure 2.2 in
%   "Computational Methods for Image Reconstruction" - Chung and Ruthotto
%       NMR Biomedicine Special Issue: QSM, 2016
%
% Chung & Ruthotto (2016)

rng(0)
%% get phantom
load simu_RDF

omega = zeros(1,6);
omega(1:2:end) = -matrix_size./2;
omega(2:2:end) = matrix_size./2;
matrix_size_pad = [256 256 128];

%% Get the eigen-values of the A matrix
Dqsm       = dipole_kernel(matrix_size_pad, voxel_size, B0_dir);
Dqsm_shift = fftshift(Dqsm);

% Investigate for Gaussian PSF
[PSF, center] = Gauss( [2 2 2] , matrix_size);
Ddb       = abs(fftn(circshift(PSF,1-center))); % Singular values
Ddb_shift = fftshift(Ddb);

%% Look at the true and observed images in the frequency domain
xtrue = iMag;
xtrue(Mask==1) = xtrue(Mask==1) - mean(xtrue(Mask==1)); 
xtrue(~Mask) = 0;

%% QSM forward operator and data
wdata = double(dataterm_mask(1,N_std,Mask)); 
Aqsm  = @(x) reshape(unpad(real(ifftn(Dqsm.*fftn(pad(reshape(x,matrix_size),matrix_size_pad)))),matrix_size),[],1);
bqsmTrue = reshape(Aqsm(xtrue),matrix_size);
bns  = .05*(randn(size(bqsmTrue))./wdata)*max(bqsmTrue(:));
bns(isnan(bns)) = 0;
bns(isinf(bns)) = 0;
bqsm = bqsmTrue + bns;
bqsm = bqsm + .3*randn(size(bqsm)).*(1-Mask);

% Deblurring operator and data
Adb  = @(x) reshape(real(ifftn(Ddb.*fftn(reshape(x,matrix_size)))),[],1);
bdbTrue  = reshape(Adb(xtrue),matrix_size);
bdb = bdbTrue + 0.05*randn(size(bdbTrue))*max(bdbTrue(:));


fig = figure(22); clf;
fig.Name = 'Fig 2.2, Chung and Ruthotto 2016: Test Data';

subplot(2,3,1)
% viewOrthoSlices2(iMag,omega,matrix_size);
imagesc(xtrue(:,:,50));
xlabel('x')
ylabel('y')
colorbar
title('ground truth')

subplot(2,3,2)
imagesc(bdb(:,:,50));
xlabel('k_x')
ylabel('k_y')
colorbar
title('simulated deblurring data')

subplot(2,3,3)
imagesc(bqsm(:,:,50));
xlabel('x')
ylabel('y')
colorbar
title('simulated QSM data')

subplot(2,3,4)
viewOrthoSlices(xtrue,omega,matrix_size);
xlabel('x')
ylabel('y')
zlabel('z')
colorbar
axis tight
title('ground truth')

subplot(2,3,5)
viewOrthoSlices(bdb,omega,matrix_size);
xlabel('x')
ylabel('y')
zlabel('z')
colorbar
axis tight
title('simulated deblurring data')

subplot(2,3,6)
viewOrthoSlices(bqsm,omega,matrix_size);
xlabel('x')
ylabel('y')
zlabel('z')
colorbar
axis tight
title('simulated QSM data')
colormap gray
