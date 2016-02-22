%% Fig51_L1.m
%  Computes the L1 reconstructions for QSM.
%  
%   This script generates L1 results shown in Figure 5.1 in
%   "Computational Methods for Image Reconstruction" - Chung and Ruthotto
%       NMR Biomedicine Special Issue: QSM, 2016
%
% Chung & Ruthotto (2016)
clear
Fig22_getTestData

delta_TE = 1;
CF = 1e6/(2*pi);
RDF = reshape(bqsm,matrix_size);
save ChungRuthottoQSMdata iMag iFreq N_std iFreq_raw RDF voxel_size matrix_size B0_dir Mask delta_TE CF xtrue
%% run SplitBregman
lambda = 40;
alpha  = 1/(2*lambda); % we want to solve alpha*TV(c) + 0.5*\|M(A*x - b)\|_2^2
[qsmSB,hisSB]     = splitBregman(1e-6,'lambda',alpha,'merit',0,'filename','ChungRuthottoQSMdata');
%% run MEDI
[qsmMEDI,hisMEDI1,hisMEDI2] = MEDI_linear('lambda',lambda,'merit',1,'filename','ChungRuthottoQSMdata');

%% show results
fig = figure(513);
fig.Name = 'Fig. 5.1. Chung and Ruthotto (2016): L1 Reconstruction';

subplot(2,3,1);
viewOrthoSlices(xtrue,omega,matrix_size);
xlabel('x')
ylabel('y')
zlabel('z')
CAX = caxis;
colorbar
axis tight

subplot(2,3,2);
viewOrthoSlices(qsmMEDI,omega,matrix_size);
xlabel('x')
ylabel('y')
zlabel('z')
colorbar
axis tight
caxis(CAX)
title('MEDI')

subplot(2,3,5);
viewOrthoSlices(abs(qsmMEDI-Mask.*xtrue),omega,matrix_size);
xlabel('x')
ylabel('y')
zlabel('z')
colorbar
axis tight
colormap gray
title(sprintf('difference, RE=%1.3f', norm(Mask(:).*(qsmMEDI(:)-xtrue(:)))/norm(Mask(:).*xtrue(:))))

subplot(2,3,3);
viewOrthoSlices(qsmSB,omega,matrix_size);
xlabel('x')
ylabel('y')
zlabel('z')
colorbar
axis tight
caxis(CAX)
title('QSM Split Bregman reconstruction')

subplot(2,3,6);
viewOrthoSlices(abs(qsmSB-Mask.*xtrue),omega,matrix_size);
xlabel('x')
ylabel('y')
zlabel('z')
colorbar
axis tight
colormap gray
title(sprintf('difference, RE=%1.3f', norm(Mask(:).*(qsmSB(:)-xtrue(:)))/norm(Mask(:).*xtrue(:))))


