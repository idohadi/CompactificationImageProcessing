function b = estimbisp(im_shc_samples, L, CGs, U, sigma)
% Estimates the spherical bispectrum of images.
% 
% Input arguments
% im_shc_samples    -   array; spherical harmonics coefficients of
%                       images. Every column is the SHC of a sampled image.
% L                 -   the band-limit of the SHCs
% U                 -   Matrix containing the debiasing coefficients
%                           U = Y_(t') * P
%                       where
%                           Y_(t')      is  evalYt(td
%                           P           is
% sigma             -   the variance of the Gaussian noise
% 
% Output arguments
%   b               -   estimator of of the spherical bispectrum
% 
% NOTE:
%   This assumes the noise is in the image space.
% 

% Compuate coefficients
UUstar = U*U';
UUT = flip(U*U.', 2);
[~, N] = size(im_shc_samples);

[b_size, ~] = size(bisp_vec(im_shc_samples(:, 1), L, CGs)); 

b = zeros(b_size, 1);
parfor j=1:N
    I_sample = im_shc_samples(:, j);
    b_im = bisp_vec(I_sample, L, CGs); 
    
    K1 = zeros(b_size, 1);
    K2 = zeros(b_size, 1);
    K3 = zeros(b_size, 1);
    
    c = 0;
    for l1=0:L
        for l2=0:l1
            for l=abs(l1-l2):l1+l2
                c = c + 1;

%               For some values of l, it exceeds the band-limit L; 
%               handling that...
                if l<=L
                    for m=-l:l
                        m11 = max(m-l2, -l1);
                        m1n = min(m+l2, l1);
                        
%                       Clebsch-Gordan ceofficients; format:
%                           cg(k) = <l1 m1(k) l2 (m-m1(k)) | l m>
%                       where 
%                           m1(k) = m11 + k - 1 
%                       for 
%                           k = 1, 2, ..., n
%                       and n is  
%                           n = m1n - m11 + 1 .
                        cg = CGs{c}{m+l+1};
                        cg = cg(:);
                        
%                       Calculate K1
                        shcIl1 = I_sample((l1^2 + 1 + m11 + l1):(l1^2 + 1 + m1n + l1));
                        shcIl1 = shcIl1(:);
                        
                        Ucoeffs = flip(UUstar(l^2 + 1 + m + l, (l2^2 + 1 + m - m1n + l2):(l2^2 + 1 + m - m11 + l2)));
                        Ucoeffs = Ucoeffs(:);
                        
                        K1(c) = K1(c) + sum(cg .* conj(shcIl1) .* Ucoeffs);
                        
%                       Calculate K2
                        shcIl2 = flip(I_sample((l2^2 + 1 + m - m1n + l2):(l2^2 + 1 + m - m11 + l2)));
                        shcIl2 = shcIl2(:);
                        
                        Ucoeffs = UUstar(l^2 + 1 + m + l, (l1^2 + 1 + m11 + l1):(l1^2 + 1 + m1n + l1));
                        Ucoeffs = Ucoeffs(:);
                        
                        K2(c) = K2(c) + sum(cg .* conj(shcIl2) .* Ucoeffs);
                        
%                       Calculate K3
                        shcIlm = I_sample(l^2 + 1 + m + l);
                        
                        rowi = l1^2 + 1 + m11 + l1;
                        coli = (L+1)^2 - (l2^2 + 1 + m - m11 + l2) + 1;
                        d = coli - rowi;
                        Ucoeffs = diag(UUT, d);
                        if d >=0
                            Ucoeffs = Ucoeffs(rowi:(l1^2 + 1 + m1n + l1));
                        else
                            Ucoeffs = Ucoeffs(coli:((L+1)^2 - (l2^2 + 1 + m - m1n + l2) + 1));
                        end
                        Ucoeffs = Ucoeffs(:);
                        
                        K3(c) = K3(c) + shcIlm * ( sum(cg .* conj(Ucoeffs)) );
                    end
                end
                
            end
        end
    end
    
    b = b + (b_im - sigma*( K1 + K2 + K3));
    
end

b = b/N;
