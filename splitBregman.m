% =========================================================================
% function  [x, his] = splitBregman(droptol,varargin)
%
% Last changed: Lars Ruthotto 2016/02/12
%
% Split Bregman method for QSM problem
%
% Input:
%  droptol   - tolerance for singular values of A (to remove rank
%              deficiency)
%  varargin  - usage as in MEDI. At least varargin={'filename',myFile.mat}
%
% Output:
%  x         - QSM reconstruction
%  his       - convergence history
% =========================================================================

function [x, his] = splitBregman(droptol,varargin)

[alpha iFreq RDF N_std iMag Mask matrix_size matrix_size0 voxel_size delta_TE CF B0_dir merit smv radius data_weighting gradient_weighting Debug_Mode] = parse_QSM_input(varargin{:});

%%%%%%%%%%%%%%% weights definition %%%%%%%%%%%%%%
lsqr_max_iter = 30;
lsqr_tol = 0.01;
max_iter = 10;
tol_norm_ratio = 0.001;
data_weighting_mode = data_weighting;
gradient_weighting_mode = gradient_weighting;
omega = zeros(1,6);
omega(2:2:end) = voxel_size.*matrix_size;
matrix_size_pad = 2.^(nextpow2(matrix_size));


L    = getGradientShortDiff(matrix_size,voxel_size);
% L    = getGradientCentered(matrix_size,voxel_size);
grad = @(x,h) L*double(x(:));
mkvec = @(x) x(:);
mkmat = @(v) reshape(v,matrix_size);
iter = 0;

if (smv)
    S = SMV_kernel(matrix_size, voxel_size,radius);
    D=S.*dipole_kernel(matrix_size, voxel_size, B0_dir);
    Mask = SMV(Mask, matrix_size,voxel_size, radius)>0.999;
    RDF = iFreq - SMV(iFreq, matrix_size, voxel_size, radius);
    RDF = RDF.*Mask;
    N_std = SMV(N_std, matrix_size, voxel_size, radius);
else
    D=dipole_kernel(matrix_size_pad, voxel_size, B0_dir);
end

w = double(mkvec(dataterm_mask(data_weighting_mode, N_std, Mask)));
RDF = double(mkvec(RDF));
wG = gradient_mask(gradient_weighting_mode, iMag, Mask, grad, voxel_size,0.9);
D(abs(D)<droptol) = 0;
F = @(x) mkvec(unpad(real(ifftn(D.*fftn(pad(mkmat(x),matrix_size_pad)))),matrix_size));
% F  = @(x) mkvec(real(ifftn(D.*fftn(mkmat(x)))));

Q = speye(prod(matrix_size));
Q = Q(:,Mask(:)==1);
wLQ = sdiag(wG)*L*Q;
% impose Neuman boundary condition
tt  = sum(L*double(Mask(:)),2);
wLQ(abs(tt)>0.01,:) = [];
% remove zero rows from weighted gradient
tt  = sum(abs(wLQ),2);
wLQ(tt==0,:) = [];

x = zeros(size(Q,2),1);
z = zeros(size(wLQ,1),1);
b = zeros(size(wLQ,1),1);
mu = 10*sqrt(alpha);
shrink = @(x,alpha) sign(x).*max(abs(x)-alpha, 0);

A  = @(x) [w.*F(Q*x); sqrt(mu)*(wLQ*x)];
At = @(x) Q'*F(w.*x(1:size(Q,1))) + sqrt(mu)*(wLQ'*x(size(Q,1)+1:end));
af = @(x,f) afun(A,At,x,f); 
    
% save some output
his = zeros(max_iter+1,8);
DOld = 0.5*norm(w.*(F(Q*x)- RDF),2)^2;
SOld = sum(abs(wLQ*x));
JOld = DOld+alpha*SOld;
his(1,:) = [JOld DOld 0 SOld 0 0 0 0];

fprintf('------------------ QSM-TV with SplitBregman --------------------------\n')
fprintf('iter\tJc\t\tDc\t\tDc-Dold\t\tSc\t\t|x-xOld|\trelres_cg\tcg_iter\ttime\n');
fprintf('%2d\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%2d\t%1.1f\n',iter,his(1,:));

xOld = x;
ndx = Inf;
while (ndx>tol_norm_ratio)&&(iter<max_iter)
tic
    iter=iter+1;
    % update x by solving least squares problem
    rhs = [w.*RDF; sqrt(mu)*(z- b)];
    [x,~,RELRES,ITER,~] = lsqr(af,rhs,lsqr_tol,lsqr_max_iter,[],[],x);
    
    % update z by using soft thresholding
    z = shrink(wLQ*x + b, alpha/mu);
    
    % update b
    b = b + (wLQ*x - z);

    iterTime = toc;
    
    % save statistics
    wres=w.*(F(Q*x)- RDF);
    Dc  = 0.5*norm(wres,2)^2;
    Sc  = sum(abs(wLQ*x));
    Jc  = Sc+alpha*Sc;
    ndx = norm(x-xOld)/norm(xOld);
    
    his(iter,:) = [Jc Dc Dc-DOld Sc ndx RELRES ITER iterTime];
    fprintf('%2d\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%2d\t%1.1f\n',iter,his(iter,:));
    % store old values
    DOld = Dc; SOld = Sc; JOld = Jc; xOld = x;
   
    
    if merit
        a = wres(Mask==1);
        ma = mean(a); factor = std(a)*5;
        wres = wres-ma;
        wres = abs(wres)/factor;
        wres(wres<1) = 1;
        N_std(Mask==1) = N_std(Mask==1).*wres(Mask==1).^2;
        w(Mask==1) = w(Mask==1)./wres(Mask==1).^2;
    end
    
    fig = figure(99);clf;
    fig.Name = sprintf('Chung and Ruthotto 2016: Split Bregman, iter=%d',iter);
    viewOrthoSlices(Q*x,omega,matrix_size);
    title(sprintf('iter=%d, intensity range=[%f,%f]',iter,min(x),max(x)));
    colormap gray
    view([52 24]);
    pause(.1);
    
    
end



%convert x to ppm
x = reshape(Q*x/(2*pi*delta_TE*CF)*1e6,matrix_size);

if (matrix_size0)
    x = x(1:matrix_size0(1), 1:matrix_size0(2), 1:matrix_size0(3));
    iMag = iMag(1:matrix_size0(1), 1:matrix_size0(2), 1:matrix_size0(3));
    RDF = RDF(1:matrix_size0(1), 1:matrix_size0(2), 1:matrix_size0(3));
    Mask = Mask(1:matrix_size0(1), 1:matrix_size0(2), 1:matrix_size0(3));
    matrix_size = matrix_size0;
end

end





              