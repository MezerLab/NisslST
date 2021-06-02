function [theta_mean, anisotropy_mean, coherence, theta_vec, thresh, mean_val] = st_get_subimg_structure_tensor_statistics(imFile,rho,sigma,sample_near_cells)
% Based on the code of this paper: https://www.math.univ-toulouse.fr/~weiss/Publis/Journals/2015/Structure_Tensor_Cell_Organization_Zhang_Weiss_2015.pdf
% Code downloaded from: https://www.math.univ-toulouse.fr/~weiss/PageCodes.html
% Path to code: /ems/elsc-labs/mezer-a/roey.schurr/CodeMatlab/functions_from_the_web/StructureTensor_toolbox
%
% Consider adding also the mean anisotropy, to plot with HSV
%
% coherence is another measure for anisotropy. It's calculated as
% follows: for each point in the image, calculate the main (structure
% tensor) orientation. Now, normalize each such vector (actually, I think
% it's already normalized). Next, calculate the sum of all vectors (this is
% like a phase consistency measure). Finally, take the norm of that vector.
% If many vectors point in opposite directions, we will get a small vector.
%
% rho=23 for AHB data: 23*0.645 =~ 15 mu
% rho=15 for abatlas data (15 mu): 15*0.97 =~ 15 mu

binarize_and_fill_holes = false; % Binarize the image using Otshu threshold before calculating structure tensor, to avoid any effect of background staining

rand('seed',0); randn('state',0)
test_plots = false; % Make some plots for testing
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isstr(imFile) % Usually we actually give as input a single small image
    im = imread(imFile);
else
    im = imFile;
end
if size(im,3)>1
    im = rgb2gray(im);
end
im = 255-im;

%% Get mean value in image
mean_val = mean(im(:));

%% Get threshold for binarizing sub-image
% This is later used to construct the brain mask
bins_centers = [0:254]+0.5;
counts = hist(im(:),bins_centers);
[thresh,~] = otsuthresh(counts);

%% Fill holes
if binarize_and_fill_holes
    im_bw = im2bw(im,thresh);
    im = imfill(im_bw,'holes');
end

%% Create sampling mask
if sample_near_cells
% This is the old version, only very close to cells. The new version is
% close up to the std of the Gaussian kernel (rho) used for spatial
% regularization
    mask = im2bw(im,thresh);
    mask = imfill(mask,'holes');
%     mask = 1-mask;
    mask = imdilate(mask,strel('disk',rho));
else
    mask = ones(size(im));
end
mask = logical(mask);

%% Compute structral tensor 
EigInfo = coherence_orientation_with_sigma(double(im),rho,sigma); 

if test_plots
    para.Step  =20; %%% intensity of orientation
    para.scl   =10; %%% length of orientation
    ConvInfo.imconv = ones(size(im));
    DisplayImage(im ,EigInfo,ConvInfo,para);
end

%% Calculate coherence of the 1st eigenvector, as a measure of anisotropy over the entire sub-image
% Using the coherence of one of the eigenvectors (they're orthogonal, so it doesn't matter which we take)

% Make all eigen vectors point in [0 180]. This is important for the
% orientation coherence calculation. Fortunately, it seems like it is
% already given that way, so the conditions below are never met

% &&& I USED W1 INSTEAD OF W2, SO I NEED TO ADD 90 DEGREES. OR SUBTRACT?
eigvec_weighted = EigInfo.w2; % w1 is in [0 180], w2 is in [-180,-90] U [90 180]. Actually, right now it's not weignted

% Mask (to near-cells, or to all the image, depending on the mask created above)
eigvec_weighted = apply_mask(eigvec_weighted,mask);


if any(~isreal(eigvec_weighted))
    theta_mean = nan;
    anisotropy_mean = nan;
    coherence = nan;
    theta_vec = nan;
    return
end

% Calculate coherence
w2_norm = sqrt(eigvec_weighted(:,1).^2+eigvec_weighted(:,2).^2);
w2 = eigvec_weighted./repmat(w2_norm,[1,2]); % each vector normalized
w2_sum = sum(w2,1);
coherence = norm(w2_sum)/(size(w2,1));

% Get theta_vec (all thetas, one per pixel in the mask)
[theta_vec,~] = cart2pol(eigvec_weighted(:,1),eigvec_weighted(:,2)); % result is in [-pi pi]
theta_vec = theta_vec*180/pi; % Now everything is in [-pi,pi]
theta_vec(theta_vec<0) = theta_vec(theta_vec<0)+180; % Now everything will be in [0,180], which is good for the colormap

%% Calculate the mean eigenvector, to get the angle
% eigvec_mean = mean(eigvec_weighted);
% Or try like this:
eigvec_mean = sum(eigvec_weighted);
eigvec_mean = eigvec_mean./norm(eigvec_mean);
[theta_mean,~] = cart2pol(eigvec_mean(1),eigvec_mean(2)); % Returns [-pi,pi]
theta_mean = theta_mean*180/pi;
if theta_mean<0
    theta_mean = theta_mean+180; % theta will be in [0,180])
end
theta_mean = round(theta_mean); % A score over the 180 colormap values

%% Calculate the mean anisotropy (aniso is local, for each tensor)
aniso = sqrt(EigInfo.nu1./EigInfo.nu2);
aniso = apply_mask(aniso,mask);
anisotropy_mean = mean(aniso);

end

%% Additional functions
function img_out = apply_mask(img,mask)
    % This will only take the values of img inside mask. Output is a
    % nx1 vector. If img has more dimensions, output is nxm.
    img_out = [];
    for dI = 1:size(img,3)
        tmp = img(:,:,dI);
        tmp = tmp(mask);
        img_out(:,dI) = tmp(:); % (nxm) X d matrix. The original img dimensions are not reduced to long vectors
    end    
end

function theta = make_between_0_180(theta)
    theta(theta<0) = theta(theta<0)+180;
    theta = round(theta);
    theta(theta==0) = 1;
    theta(theta==181) = 180;
end

