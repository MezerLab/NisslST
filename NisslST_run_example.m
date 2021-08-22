clear;clc
%% Prerequisites
% Download the MATLAB toolbox for Structure Tensor Analysis from:
% https://www.math.univ-toulouse.fr/~weiss/PageCodes.html
% and define its location here:
st_toolbox_dir = 'ADD/CODE/TO/THE/STRUCTURE_TENSOR_TOOLBOX/DIRECTORY';

%% User settings
nissl_st_dir = 'ADD/CODE/TO/THE/NisslST/DIRECTORY';
addpath(genpath(nissl_st_dir));
addpath(genpath(st_toolbox_dir));

%% Load the matplotlib colormaps
matplotlib_colormaps

%% Download a histological image of the human corpus callosum from Allen Brain Atlas
% This should download the example jpg to your default downloads folder.
im_id = '112360908';
top = 35000;
left = 5000;
w = 20000;
h = 17000;

cmd = sprintf('http://api.brain-map.org/api/v2/image_download/%s?left=%g&top=%g&width=%g&height=%g',im_id,left,top,w,h);
web(cmd)

%% Set paths
imFile = 'ADD/FULL/PATH/TO/THE/DOWNLOADED/IMAGE';

output_dir = fullfile(nissl_st_dir,'output');
if ~exist(output_dir,'dir')
    mkdir(output_dir)
end
figDir = fullfile(output_dir,'figures');
if ~exist(figDir,'dir')
    mkdir(figDir)
end


%% Set analysis flags and options
plot_flag = true;
save_plots = false;
ds_factor = 1;% Downsample factor for downsampling the image
pix_microns = 0.97;
side_microns = 200; % Tile side in microns
rho = 15; % measured in pixels (~15 microns)
sigma = 3; % Gaussian kernel for smoothing the image before calculating Nissl-ST. This parameter is typically 0, but for this Allen Brain Atlas dataset, it seems that sigma=3 yields better results.
nPeaks = 2; % Number of peak orientations to extract. In fact, we typically just take the first peak.
sample_near_cells = true; % Only extract orientations as close as 'rho' to cells (to avoid measuring empty spaces)

%% Apply Nissl-ST and and save analysis results to file
% The st_run_on_medium_img() function runs on a medium size image
% (typically less than a whole slice, which might be too large to load into
% memory, and including several thousands of small image tiles (around
% 50-200^2 microns^2). It saves the analysis output to a file (st_stats_file)
% that can later be loaded for further anaylsis and visualization.
%
% This should take less than an hour for 200 microns tiles.

tic
st_stats_file = fullfile(output_dir,sprintf('st_stats_side_%gmu_rho%g_sigma%g_npeaks%g_samplenearcells_%g_slicethresh_ds%g.mat',side_microns,rho,sigma,nPeaks,sample_near_cells,ds_factor));
if ~exist(st_stats_file,'file')
    st_run_on_medium_img(imFile,pix_microns,side_microns,rho,sigma,nPeaks,sample_near_cells,ds_factor,st_stats_file);
end
toc

%% Display a single image tile with overlaid local tensors 
im = imread(imFile);
im = rgb2gray(im);
side_pix = round(side_microns/pix_microns); % The side of each sub-image in pixels

upper_left_x = 100;
upper_left_y = 14500;
rect = [upper_left_x, upper_left_y, side_pix, side_pix];
sub_img = imcrop(im,rect);

plot_version = true;
st_get_img_structure_tensor_statistics(sub_img,side_pix,rho,sigma,nPeaks,sample_near_cells,ds_factor,plot_version)

%% Load the structure tensor statistics file
load(st_stats_file);

%% Reorient all maps
mean_val_map = reorient_img(mean_val_map);
thresh_map = reorient_img(thresh_map);
theta_mean_map = reorient_img(theta_mean_map);
aniso_map = reorient_img(aniso_map);
coherence_map = reorient_img(coherence_map);
coherence_map = coherence_map.^2;
theta_peaks_map = reorient_img(theta_peaks_map);
peaks_height_map = reorient_img(peaks_height_map);
theta_von_mises_map = reorient_img(theta_von_mises_map);
von_mises_component_fraction_map = reorient_img(von_mises_component_fraction_map);
cluster_pixel_size_map = gm_wm_info.cluster_pixel_size;
cluster_pixel_size_map = reorient_img(cluster_pixel_size_map);

%% Get initial brain mask based on the Otsu threshold per sub-image
% Note: The functions used for creating the masks typically work best when
% using entire slices.
mask = getMaskFromThreshMap(thresh_map);

% Get GM mask based on the distribution of cluster cell sizes in pixels (look at 95th prctile)
[gmMask, wmMask] = get_gm_wm_masks(cluster_pixel_size_map,mask);

% Get a grayscale image with GM/WM contrast
grayim = squeeze(cluster_pixel_size_map(:,:,1));
grayim(abs(grayim)>400) = 0;
grayim(~mask) = 0;

% Convert angle maps to rgb images
theta_peaks_rgb = theta_to_rgb(theta_peaks_map,mask);

%% Get limits for the grayscale background image
grayim = mean_val_map;
vals = grayim(mask);
clims = prctile(vals,[1,99]);
grayim(grayim>clims(2)) = clims(2);
grayim(grayim<clims(1)) = clims(1);
grayim = grayim./mean(grayim(grayim>0));
grayim = grayim.*mask;
grayim = grayim./max(grayim(:));
imshow(grayim)
vals = grayim(mask);
clims = prctile(vals,[1,99]);

%% Plot peak RGB orientation
vI = 1; % Plot the first peak
figure('Color','k')
rgb_shaded = theta_peaks_rgb{vI}.*repmat(coherence_map,[1,1,3]);
imshow_rgb_on_gray(rgb_shaded,grayim,logical(wmMask),clims) % Change clims_all to clims if not working on downsampled images
axis equal
if save_plots
    set(gcf,'position',get(0,'screensize'))
    export_fig(gcf,fullfile(figDir,sprintf('%s_ornttn_rgb_shdd_by_coherence_side_%g_rho_%g_sigma_%g_sample_near_cells_%g_crtx_peak_%g_noshading.png', side_microns, rho, sigma, sample_near_cells, vI)),'-dpng','-r300');
end

%% Plot coherence map
tmp = coherence_map;
tmp(~mask) = 0;
coherence_plot = ind2rgb(round(tmp*255),infernodata);
figure('Color','k')
imshow_rgb_on_gray(coherence_plot,grayim,logical(wmMask),clims)
axis equal
if save_plots
    set(gcf,'position',get(0,'screensize'))
    export_fig(gcf,fullfile(figDir,sprintf('%s_cohrence_side_%g_rho_%g_sigma_%g_sample_near_cells_%g_crtx.png',sliceName, side_microns, rho, sigma, sample_near_cells)),'-dpng','-r300');
end
