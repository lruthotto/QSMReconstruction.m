function [PSF, center] = Gauss( sigma, N)
%
%     [PSF, center] = Gauss( sigma, N)
%
% This function constructs the Gaussian blur PSF, and the center of the PSF
% for 1, 2, and 3D blurs.
%
%   Inputs:
%    sigma - parameters defining the spread of the blur
%        N - size of the PSF, determines the dimension of the PSF
%
%    Output:
%      PSF - normalized Gaussian point spread function
%            For 1D: exp(- .5 x^2 / sigma^2)
%            For 2D: 
%               exp(-x^2/2sigma^2 - y^2/2sigma^2), if sigma is a number  
%               exp(-x^2/2sigma(1)^2 - y^2/2sigma(2)^2), if sigma is 
%                                                           length 2 vector
%               elliptical Gaussian is sigma is length 3 vector
%            For 3D:
%               exp(-x^2/2sigma(1)^2 - y^2/2sigma(2)^2- z^2/2sigma(3)^2), 
%                   if sigma is length 3 vector
%  center  - center of the PSF
%
%  J. Chung 8/2015
%

dim = length(N);

switch dim
    case 1 % 1D PSF
        x = -fix(N/2):ceil(N/2)-1;
        PSF = exp( -.5*(x.^2)/(sigma^2));
        PSF = PSF/sum(PSF(:));
        center = N/2+1;
    
    case 2 % 2D PSF
        m = N(1); n = N(2);
        
        % Set up grid points to evaluate the Gaussian function.
        x = -fix(n/2):ceil(n/2)-1;
        y = -fix(m/2):ceil(m/2)-1;
        [X,Y] = meshgrid(x,y);
        
        % Compute the Gaussian PSF
        if length(sigma) == 1
            PSF = exp( -(X.^2)/(2*sigma^2) - (Y.^2)/(2*sigma^2) );
            PSFsum = sum(PSF(:));
            
        elseif length(sigma) == 2
            s1 = sigma(1); s2 = sigma(2);
            PSF = exp( -(X.^2)/(2*s1^2) - (Y.^2)/(2*s2^2) );
            PSFsum = sum(PSF(:));
            
        elseif length(sigma)==3 % Elliptical Gaussian function
            s1 = sigma(1); s2 = sigma(2); s3 = sigma(3);
            num = -((X.^2)*(s1^2) + (Y.^2)*(s2^2) - 2*(s3^2)*(X.*Y));
            den = 2*(s1^2 * s2^2 - s3^4);
            PSF = exp( num / den );
            PSFsum = sum(PSF(:));
            
        end
        PSF = PSF/PSFsum; % normalize the PSf
        center = [m/2+1, n/2+1];
    
    case 3 % 3D PSf
        % N is a vector with [m,n,k]
        m = N(1); n = N(2); k = N(3);
        
        % Set up grid points to evaluate the Gaussian function.
        x = -fix(n/2):ceil(n/2)-1;
        y = -fix(m/2):ceil(m/2)-1;
        z = -fix(k/2):ceil(k/2)-1;
        [X,Y,Z] = ndgrid(x,y,z);
        
        % Compute the Gaussian PSF
        s1 = sigma(1); s2 = sigma(2); s3 = sigma(3);
        PSF = exp( -(X.^2)/(2*s1^2) - (Y.^2)/(2*s2^2) - (Z.^2)/(2*s3^2) );
        PSFsum = sum(PSF(:));
        
        PSF = PSF/PSFsum;
        center = [m/2+1, n/2+1, k/2+1];
        
end
