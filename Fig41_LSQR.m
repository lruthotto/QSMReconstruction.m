%% Fig41_LSQR2.m
%  Use LSQR to illustrate semi-convergence for the image deblurring and QSM
%  problems.
%  
%   This script generates Figure 4.1 in
%   "Computational Methods for Image Reconstruction" - Chung and Ruthotto
%       NMR Biomedicine Special Issue: QSM, 2016    
%
% Chung & Ruthotto (2016)

clear
Fig22_getTestData

maxIter = 20; % number of LSQR iterations

%% Solve the QSM problem using LSQR
wdata = wdata(:); % weights for data term

% get QSM forward operator
threshQSM = 1e-4;
Dt  = Dqsm;
Dt(abs(Dt)<threshQSM) = threshQSM;
data = unpad(real(ifftn(Dt.*fftn(pad(bqsm,matrix_size_pad)))),matrix_size);
Aqsm   = @(x) reshape(unpad(real(ifftn(Dt.*fftn(pad(reshape(x,matrix_size),matrix_size_pad)))),matrix_size),[],1);
wmask = Mask(:);


P = speye(prod(matrix_size));
P = P(:,Mask(:)==1);

A = @(x) wdata.*Aqsm(P*x); % forward
At = @(x) P'*Aqsm(wdata.*x); % adjoint
Afun3 = @(v,flag) afun(A,At,v,flag);

%% QSM: Run LSQR for noise free data
input = HyBRset('Iter',maxIter, 'InSolv', 'none', 'x_true', P'*xtrue(:));
[~, output_lsqr] = HyBR(Afun3, wdata.*double(bqsmTrue(:)), [], input, 1.0);
errQSM1(:,1) = output_lsqr.Enrm;
errQSM1(:,2) = output_lsqr.Rnrm;

%%
fig = figure(41);clf;
fig.Name = 'Fig 4.1, Chung and Ruthotto 2016: Semiconvergence';
subplot(2,2,2)
[mE,iE] = min(errQSM1(:,1)); 
plotyy(1:21,errQSM1(:,1),1:maxIter+1,errQSM1(:,2));
hold on
plot(iE,mE,'or');
legend('||xtrue-x_k||','x_{opt}','misfit')
title('QSM, noise free data')

%% QSM:  Run LSQR for noisy data
input = HyBRset('Iter',maxIter, 'InSolv', 'none', 'x_true', P'*xtrue(:));
[~, output_lsqr] = HyBR(Afun3, wdata.*double(bqsm(:)), [], input, 1.0);
errQSM2(:,1) = output_lsqr.Enrm;
errQSM2(:,2) = output_lsqr.Rnrm;

%%
subplot(2,2,4)
[mE,iE] = min(errQSM2(:,1)); 
ax = plotyy(1:21,errQSM2(:,1),1:maxIter+1,errQSM2(:,2));
hold(ax(1),'on')
hold(ax(2),'on')
plot(ax(1),iE,mE,'or');
nl = norm(sqrt(wdata(:)).*double(bqsm(:)-bqsmTrue(:)));
plot(ax(2),1:maxIter+1,nl*ones(maxIter+1,1),'--k');

legend('||xtrue-x_k||','x_{opt}','misfit','noise level')
title('QSM, noisy data')

%% Deblurring: Run LSQR for noise free data
A = @(v,flag) afun(Adb,Adb,v,flag);

input = HyBRset('Iter',maxIter, 'InSolv', 'none', 'x_true', xtrue(:));
[~, output_lsqr] = HyBR(A, bdbTrue(:), [], input, wmask);
errdb1(:,1) = output_lsqr.Enrm;
errdb1(:,2) = output_lsqr.Rnrm;

%%
subplot(2,2,1)
[mE,iE] = min(errdb1(:,1)); 
plotyy(1:21,errdb1(:,1),1:maxIter+1,errdb1(:,2));
hold on
plot(iE,mE,'or');
% legend('||xtrue-x_k||','x_{opt}','misfit')
title('Deblurring, noise free data')

%% Deblurring: Run LSQR for noisy data

input = HyBRset('Iter',maxIter, 'InSolv', 'none', 'x_true', xtrue(:));
[~, output_lsqr] = HyBR(A, bdb(:), [], input, wmask);
errdb2(:,1) = output_lsqr.Enrm;
errdb2(:,2) = output_lsqr.Rnrm;

%%
subplot(2,2,3)
[mE,iE] = min(errdb2(:,1)); 
ax =plotyy(1:21,errdb2(:,1),1:maxIter+1,errdb2(:,2));
hold(ax(1),'on')
hold(ax(2),'on')
plot(ax(1),iE,mE,'or');
nl = norm(bdb(:)-bdbTrue(:)); % noise level
plot(ax(2),1:maxIter+1,nl*ones(maxIter+1,1),'--k');
% legend('||xtrue-x_k||','x_{opt}','misfit')
title('Deblurring, noisy data')

