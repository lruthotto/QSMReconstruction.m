%% Fig51_TSVD.m
%  Computes the TSVD reconstruction for QSM.
%  
%   This script generates TSVD results shown in Figure 5.1 in
%   "Computational Methods for Image Reconstruction" - Chung and Ruthotto
%       NMR Biomedicine Special Issue: QSM, 2016
%
% Chung & Ruthotto (2016)

%% get test data
Fig22_getTestData

% remove background where there is no reliable data
bqsm = bqsm.*Mask;
bqsmTrue = bqsmTrue.*Mask;

%% Compute TSVD solution
discr = norm(bqsmTrue(:)-bqsm(:)); % norm of noise
alpha  = 2.3e-1;                    % tolerance for singular value
D_tsvd = Dqsm(:);
ids = abs(D_tsvd) < alpha;
D_tsvd(ids) = 0;
D_tsvd(~ids) = 1./D_tsvd(~ids);
D_tsvd = reshape(D_tsvd,matrix_size_pad);
bh = fftn(pad(reshape(bqsm,matrix_size),matrix_size_pad));
x_tsvd = unpad(real(ifftn(bh.*D_tsvd)),matrix_size);

b_tsvd = Aqsm(x_tsvd);
err    = norm(b_tsvd(:)-bqsm(:));
fprintf('err=%1.2e vs discr=%1.2e\n',err,discr);
%%
fig = figure(511); clf;
fig.Name = 'Fig 5.1, Chung and Ruthotto (2016): TSVD Reconstruction for QSM';
subplot(1,2,1)

viewOrthoSlices2D(flip(reshape(xtrue,matrix_size),1),[0 1 0 1 0 1],matrix_size);
title('true image');
cax = caxis;
colormap gray

subplot(1,2,2);
viewOrthoSlices2D(flip(reshape(x_tsvd,matrix_size),1),[0 1 0 1 0 1],matrix_size);
title(sprintf('TSVD reconstruction, alpha=%1.3f, err=%1.3e, discrepancy=%1.3e', alpha,err,discr));
colormap gray
caxis(cax)
