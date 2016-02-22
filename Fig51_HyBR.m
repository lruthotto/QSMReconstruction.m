%% Fig51_HyBR.m
%  Perform QSM reconstruction using hybrid iterative regularization
%  
%   This script generates Figure 5.1 in
%   "Computational Methods for Image Reconstruction" - Chung and Ruthotto
%       NMR Biomedicine Special Issue: QSM, 2016   
%
% Chung & Ruthotto (2016)

Fig22_getTestData;

%% 

% get projector from brain mask to whole image
P = speye(prod(matrix_size));
P = P(:,Mask==1);

% get ingredients of forward problem
D  = dipole_kernel(matrix_size_pad, voxel_size, B0_dir);
A  = @(x) reshape(unpad(real(ifftn(D.*fftn(pad(reshape(P*x,matrix_size),matrix_size_pad)))),matrix_size),[],1);

% remove rank deficiency
tol = 1e-5;
Dt  = D;
Dt(abs(Dt)<tol) = 0;
At = @(x) P'*reshape(unpad(real(ifftn(Dt.*fftn(pad(reshape(x,matrix_size),matrix_size_pad)))),matrix_size),[],1);
bqsm = reshape(unpad(real(ifftn((abs(D)>tol).*fftn(pad(reshape(bqsm,matrix_size),matrix_size_pad)))),matrix_size),[],1);

maxIter = 20;

%%  solve least squares problem with masking + data weighting
wmask = Mask(:);
wdata = wdata(:);
Af3 = @(x) wdata.*A(x);
At3 = @(x) At(wdata.*x);
Afun3 = @(v,flag) afun(Af3,At3,v,flag);

% LSQR
input = HyBRset('Iter',maxIter, 'InSolv', 'none', 'x_true', P'*xtrue(:));
[x_lsqrt, output_lsqr] = HyBR(Afun3, wdata.*bqsm(:), [], input, 1.0);
x_lsqr = P*x_lsqrt;
Enrm_lsqr = output_lsqr.Enrm;

% HyBR with optimal regularization parameter
input = HyBRset('Iter',maxIter,'RegPar', 'optimal', 'x_true', P'*xtrue(:));
[x_optt, output_opt] = HyBR(Afun3, wdata.*bqsm(:),[], input, 1.0);
x_opt = P*x_optt;
Enrm_opt = output_opt.Enrm;

% HyBR with GCV parameter
input = HyBRset('Iter',maxIter, 'RegPar', 'gcv', 'x_true', P'*xtrue(:));
[x_gcvt, output_gcv] = HyBR(Afun3, wdata.*bqsm(:),[], input, 1.0);
x_gcv = P*x_gcvt;
Enrm_gcv = output_gcv.Enrm;

% HyBR with DP parameter
nlvl = norm(wdata.*bqsm(:)-Afun3(P'*xtrue(:),'notransp'));
input = HyBRset('Iter',maxIter,'RegPar', 'DP', 'nLevel', nlvl, 'x_true', P'*xtrue(:));
[x_dpt, output_dp] = HyBR(Afun3, wdata.*bqsm(:), [],input, 1.0);
x_dp = P*x_dpt;
Enrm_dp = output_dp.Enrm;

err = [Enrm_lsqr(:), Enrm_opt(:), Enrm_gcv(:), Enrm_dp(:)];

%% Display Results
fontSize = 30;
fig = figure(512);
fig.Name = 'Fig 5.1, Chung and Ruthotto 2016: HyBR';
plot(err(:,1),'LineWidth',2), hold on
plot(err(:,2), 'r--','LineWidth',2),
plot(err(:,3), 'k:','LineWidth',2),
plot(err(:,4), 'm-.','LineWidth',2)
plot(20, err(20,3), 'k*','MarkerSize',10,'LineWidth',2)
xlabel('iteration, k')
ylabel('relative error')
legend('LSQR','HyBR OPT', 'HyBR GCV', 'HyBR DP')
axis([1,maxIter,.55,.8])
set(gca,'FontSize',fontSize)

% Display reconstructions
fig = figure(2); clf;
subplot(1,3,1),
imagesc(xtrue(:,:,50)); title('True Image')
CAXtrue = caxis;
xlabel('x')
ylabel('y')
axis image
colormap gray;

% Solution at 20 iterations of LSQR
subplot(1,3,2),
x_lsqr = reshape(x_lsqr,matrix_size);
imagesc(x_lsqr(:,:,50));
caxis(CAXtrue);
colormap gray
axis image
title('LSQR')

% Solution at 20 iterations of HyBR GCV
subplot(1,3,3),
x_gcv = reshape(x_gcv,matrix_size);
imagesc(x_gcv(:,:,50));
caxis(CAXtrue);
colormap gray
axis image
title('HyBR GCV')

% Error images
eimage1 = abs(x_lsqr(:,:,50) - xtrue(:,:,50));
eimage2 = abs(x_gcv(:,:,50) - xtrue(:,:,50));
figure(3)
subplot(1,2,1), imshow(eimage1,[]), colormap(flipud(colormap))
subplot(1,2,2), imshow(eimage2,[]), colormap(flipud(colormap))
