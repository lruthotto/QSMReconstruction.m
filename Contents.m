%% Contents.m
%   This folder contains MATLAB code to accompany the invited review paper:
% 
%   "Computational Methods for Image Reconstruction" - Chung and Ruthotto
%       NMR Biomedicine Special Issue: QSM, 2016    
%
% These codes require the MEDI package from the Cornell MRI research lab
% that can be obtained from: http://weill.cornell.edu/mri/pages/qsmreview.html
%
% Chung & Ruthotto (2016)
%% Script Files used to generate figures in the paper
%
%  Fig22_getTestData.m    Generates an image deblurring example and a QSM
%                           reconstruction example. (Fig 2.2)
%
%  Fig21_PSFs.m           Compares the point spread function for image 
%                           deblurring and QSM. (Fig 2.1)
%
%  Fig31_DiscretePicard   Computes the discrete Picard plots for image
%                           deblurring and QSM. (Fig 3.1)
%   
%  Fig41_LSQR.m           Use LSQR to illustrate semi-convergence for 
%                           image deblurring and QSM. (Fig 4.1)
%
%  Fig51_HyBR.m           Perform QSM reconstruction using hybrid iterative
%                           regularization. (Fig 5.1)
%
%  Fig52_L1               Perform QSM reconstruction using l_1
%                           regularization (MEDI and Split Bregman). 
%                           (Fig 5.2)
%
%  Fig52_TSVD            Perform TSVD reconstruction for QSM problem
%                           (Fig 5.2)
%
%% MATLAB functions used to generate figures in the paper
%
%  afun.m                Helper function for evaluating linear operator and
%                           transpose
%
%  Gauss.m               Helper function constructing a Gaussian blurring
%                           kernel
%
%  getGradientShortDiff.m  Generates short difference discretization of
%                             gradient operator
%
%  pad.m                Helper function to pad an image
%
%  sdiag.m                Helper function for generating a sparse diagonal
%                              matrix
%
%  splitBregman.m       Split Bregman implementation for QSM
%
%  unpad.m              Helper function to remove padding from image
%
%  viewOrthoSlices2D.m  Visualizes 3D volume along 2D slices
%
%  viewOrthoSlices.m    Visualizes 3D volume along 2D slices
%

