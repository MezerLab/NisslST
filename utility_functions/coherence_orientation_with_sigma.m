function EigInfo=coherence_orientation_with_sigma(im,rho,sigma) % Roey added sigma here     
%   This function generates the eigenvalues and the eigenvectors of the
%   2D structure tensor matrix.
%
%   INPUT:
%   im: the test m-by-n 2D image
%   rho: the parameter for the filter K_\rho
%
%   OUTPUT:
%   nu1,nu2: n x m matrices containing respectively the largest and the 
%   smallest eigenvalue of the structure tensor at each pixel of the image.
%   
%   w1,w2: n x m x 2 matrices containig the eigenvectors corresponding to
%   respectivly the largest and smallest eigenvalues of structure tensor.

% %  Written by Dr. Wenxing Zhang, Omar Dounia and Pierre Weiss, 
% %  ITAV, Toulouse, France, 2014.
% %  Email: pierre.armand.weiss@gmail.com; wenxing84@gmail.com.

% % This version of the function was modified by Roey Schurr, HUJI, Isrel, 2021.
% % Email: roey.schurr@mail.huji.ac.il. It includes sigma as an input argument. sigma is the width of the Gaussian kernel that smoothes the image before the structure tensor is calculated. Use sigma=0 for no smoothing.

im = im/max(im(:));
[J_rho,im_sigma]=structure_tensor(im,sigma,rho);


Trace= J_rho(:,:,1)+J_rho(:,:,4);
Det  = J_rho(:,:,1).*J_rho(:,:,4)-J_rho(:,:,2).*J_rho(:,:,3);
DETA = sqrt(Trace.^2-4*Det); %%%% Delta

nu1  = (Trace-DETA)/2;  %%% smaller eigvalue and eigvector
w1(:,:,1) = -J_rho(:,:,2);   w1(:,:,2) = J_rho(:,:,1)-nu1;

nu2 = (Trace+DETA)/2;  %%% larger eigvalue and eigvector
w2(:,:,1) = J_rho(:,:,4)-nu2;  w2(:,:,2) = -J_rho(:,:,3); 

%%%%%%%%%%%%%%%%% normalization of eigvector
Tep = sqrt(w1(:,:,1).^2+w1(:,:,2).^2);
w1  = w1./repmat(Tep,[1,1,2]);
Tep = sqrt(w2(:,:,1).^2+w2(:,:,2).^2);
w2  = w2./repmat(Tep,[1,1,2]);

clear DETA  Det  Trace  im  Tep

EigInfo.nu1 = nu1; EigInfo.w1 = w1;
EigInfo.nu2 = nu2; EigInfo.w2 = w2;
EigInfo.J_rho   = J_rho;

end
