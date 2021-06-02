clear;clc
matplotlib_colormaps;

%% User settings
nissl_st_dir = '/ems/elsc-labs/mezer-a/code/roey/NisslST'; %'ADD/CODE/TO/THE/NisslST/DIRECTORY';
addpath(genpath(nissl_st_dir));
addpath(genpath('/ems/elsc-labs/mezer-a/roey.schurr/CodeMatlab/functions_from_the_web/StructureTensor_toolbox'));

%% Set paths
imFile = fullfile(nissl_st_dir,'example_data','AHB_slice_R3_0570_CC.png');
output_dir = fullfile(nissl_st_dir,'output');
if ~exist(output_dir,'dir')
    mkdir(output_dir)
end
figDir = fullfile(output_dir,'figures');


%% Set analysis flags and options
plot_flag = true;
save_plots = false;
binary_flag = 0; % This only holds for the naming of files. You actually need to change it within the st code
ds_factor = 1;% Downsample factor for downsampling the image
pix_microns = 0.645;
side_microns = 100; % Tile side in microns
rho = 23; % measured in pixels (23 pixels of 0.645 microns = 15 microns)
sigma = 0;
nPeaks = 2; % Number of peak angles (either from findpeaks or from a von Mises mixture fit)
sample_near_cells = true;


%% Calculate and save to file
% This should take ~15 minutes for 200 microns tiles,
% nd ~30 minutes for 50 microns tiles

tic
st_stats_file = fullfile(output_dir,sprintf('ahb_st_stats_side_%gmu_rho%g_sigma%g_npeaks%g_samplenearcells_%g_slicethresh_ds%g.mat',side_microns,rho,sigma,nPeaks,sample_near_cells,ds_factor)); % &&& Change after downsmaple test
if ~exist(st_stats_file,'file')
    st_run_on_medium_img(imFile,pix_microns,side_microns,rho,sigma,nPeaks,sample_near_cells,ds_factor,st_stats_file);
end
toc

%% Load the structure tensor statistics file
load(st_stats_file);

%% Reorient all maps
mean_val_map = reorient_img(mean_val_map);
thresh_map = reorient_img(thresh_map);
theta_mean_map = reorient_img(theta_mean_map);
aniso_map = reorient_img(aniso_map);
coherence_map = reorient_img(coherence_map);
coherence_map = coherence_map.^2; % For better visualization
theta_peaks_map = reorient_img(theta_peaks_map);
peaks_height_map = reorient_img(peaks_height_map);
theta_von_mises_map = reorient_img(theta_von_mises_map);
von_mises_component_fraction_map = reorient_img(von_mises_component_fraction_map);
cluster_pixel_size_map = reorient_img(gm_wm_info.cluster_pixel_size);

%% Get mask based on the Otsu threshold per sub-image
mask = getMaskFromThreshMap(thresh_map);

% Get GM mask based on the distribution of cluster cell sizes in pixels (look at 95th prctile)
[gmMask, wmMask] = get_gm_wm_masks(cluster_pixel_size_map,mask);

% Get a grayscale image with GM/WM contrast
grayim = squeeze(cluster_pixel_size_map(:,:,1));
grayim(abs(grayim)>400) = 0;
grayim(~mask) = 0;


% Convert angle maps to rgb images
theta_mean_rgb = theta_to_rgb(theta_mean_map,mask);
theta_peaks_rgb = theta_to_rgb(theta_peaks_map,mask);
theta_von_mises_rgb = theta_to_rgb(theta_von_mises_map,mask);


%% Get limits for the grayscale background image
grayim = mean_val_map;
vals = grayim(mask);
clims = prctile(vals,[1,99]);
grayim(grayim>clims(2)) = clims(2);
grayim(grayim<clims(1)) = clims(1);
grayim = grayim./mean(grayim(grayim>0));
grayim = grayim.*mask;
imshow(grayim)
vals = grayim(mask);
clims = prctile(vals,[1,99]);
clims(2) = clims(2)*1.5; % Make image a bit darker

%% Plot the mean grayscale value
figure('Color','k')
imshow(grayim./max(grayim(:)))
% caxis(prctile(vals,[1,90]))


%% Plot peak RGB orientation
vI = 1; % Plot the first peak
figure('Color','k')
rgb_shaded = theta_peaks_rgb{vI}.*repmat(coherence_map,[1,1,3]);
imshow_rgb_on_gray(rgb_shaded,grayim./1.7,logical(mask),clims) % Change clims_all to clims if not working on downsampled images
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
imshow_rgb_on_gray(coherence_plot,grayim,logical(mask),clims)
axis equal
if save_plots
    set(gcf,'position',get(0,'screensize'))
    export_fig(gcf,fullfile(figDir,sprintf('%s_cohrence_side_%g_rho_%g_sigma_%g_sample_near_cells_%g_crtx.png',sliceName, side_microns, rho, sigma, sample_near_cells)),'-dpng','-r300');
end
