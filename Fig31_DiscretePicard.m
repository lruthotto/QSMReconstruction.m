%% Fig31_DiscretePicard.m
%  Computes the discrete Picard plots for image deblurring and QSM 
%   reconstruction.
%  
%   This script generates Figure 3.1 in
%   "Computational Methods for Image Reconstruction" - Chung and Ruthotto
%       NMR Biomedicine Special Issue: QSM, 2016    
%
% Chung & Ruthotto (2016)

%% get test data
Fig22_getTestData
bqsm = bqsm.*Mask;
bqsmTrue = bqsmTrue.*Mask;
%% Get discrete Picard plot for QSM (with noise)
Dqsm_abs = abs(Dqsm); 
Dqsm_abs = Dqsm_abs(:);
[s_qsm,idx] = sort(Dqsm_abs(:), 'descend');

% With Noise
bhat = abs(fftn(pad(reshape(bqsm,matrix_size),matrix_size_pad)));
bhat = bhat(:);
bsort_qsm = bhat(idx);

% No noise
bhat_no = abs(fftn(pad(reshape(bqsmTrue,matrix_size),matrix_size_pad)));
bhat_no = bhat_no(:);
bsort_qsm_no = bhat_no(idx);

fig = figure(1); clf;
fig.Name = 'Fig 3.1, Chung and Ruthotto 2016: Discrete Picard';
subplot(2,2,1);
skip=5000;
semilogy(1:skip:numel(bsort_qsm),s_qsm(1:skip:end),'-.k')
hold on, semilogy(1:skip:numel(bsort_qsm),bsort_qsm(1:skip:end),'ob')
semilogy(1:skip:numel(bsort_qsm),bsort_qsm(1:skip:end)./s_qsm(1:skip:end),'d', 'Color',[0 .5 0])
axis tight
% legend('|\lambda_i|','|f_i^* b|','|f_i^* b|/|\lambda_i|','Location','SouthWest')
title(sprintf('QSM, discrete Picard plot, every %d entries',skip))

%% Zoomed version with linear fit (with noise and no noise)
subplot(2,2,2);
skip=1000;
disp = 3e5;
semilogy(1:skip:disp,s_qsm(1:skip:disp),'-.k') % Singular values
hold on, semilogy(1:skip:disp,bsort_qsm(1:skip:disp),'ob') % with noise
hold on, semilogy(1:skip:disp,bsort_qsm_no(1:skip:disp),'xr') % no noise
axis tight

% Get linear fits
format long
fit = 5e4; 
p1 = polyfit(1:fit,log(s_qsm(1:fit))',1);
xvals = linspace(1,disp,100);
% yvals1 = polyval(p1,xvals);
% hold on, semilogy(xvals,exp(yvals1),'m--','LineWidth',2)

p2 = polyfit(1:fit,log(bsort_qsm(1:fit))',1);
yvals2 = polyval(p2,xvals);
hold on, semilogy(xvals,exp(yvals2),'b-','LineWidth',2) % With Noise

p2_no = polyfit(1:fit,log(bsort_qsm_no(1:fit))',1);
yvals2_no = polyval(p2_no,xvals);
semilogy(xvals,exp(yvals2_no),'r-','LineWidth',2) % No Noise
title('QSM, zoomed version with linear with (with and without noise)')
%% Get a discrete picard plot for deblurring
Ddb_abs = abs(Ddb); 
Ddb_abs = Ddb_abs(:);
[s,idx] = sort(Ddb_abs(:), 'descend');

% With Noise
bhat = abs(fftn(reshape(bdb,matrix_size)));
bhat = bhat(:);
bsort_db = bhat(idx);

% No noise
bhat_no = abs(fftn(reshape(bdbTrue,matrix_size)));
bhat_no = bhat_no(:);
bsort_db_no = bhat_no(idx);

subplot(2,2,3)
skip = 5000;
semilogy(1:skip:numel(bsort_db),s(1:skip:end),'-.k')
hold on, semilogy(1:skip:numel(bsort_db),bsort_db(1:skip:end),'ob')
semilogy(1:skip:numel(bsort_db),bsort_db(1:skip:end)./s(1:skip:end),'d','Color',[0 .5 0])
axis tight
legend('|\lambda_i|','|f_i^* b|','|f_i^* b|/|\lambda_i|','Location','SouthWest')
title(sprintf('Deblurring, discrete Picard plot, every %d entries',skip))

%% Zoomed version with linear fit

subplot(2,2,4)
skip = 1000;
disp = 3e5;
semilogy(1:skip:disp,s(1:skip:disp),'-.k')
hold on, semilogy(1:skip:disp,bsort_db(1:skip:disp),'ob') % with noise
hold on, semilogy(1:skip:disp,bsort_db_no(1:skip:disp),'xr') % no noise
axis tight

% Linear fits
format long
fit = 5e4; 
p1 = polyfit(1:fit,log(s(1:fit))',1);
xvals = linspace(1,disp,100);
% yvals1 = polyval(p1,xvals);
% hold on, semilogy(xvals,exp(yvals1),'m--','LineWidth',2)

p2 = polyfit(1:fit,log(bsort_db(1:fit))',1);
yvals2 = polyval(p2,xvals);
hold on, semilogy(xvals,exp(yvals2),'b-','LineWidth',2) % with noise

p2_no = polyfit(1:fit,log(bsort_db_no(1:fit))',1);
yvals2_no = polyval(p2_no,xvals);
hold on, semilogy(xvals,exp(yvals2_no),'r-','LineWidth',2) % no noise
title('Deblurring, zoomed version with linear with (with and without noise)')


legend('|\lambda_i|','|f_i^* b|','|f_i^* b_{true}|','Location','SouthWest')

