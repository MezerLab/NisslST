clear; clc

%% User settings
nissl_st_dir = 'ADD/CODE/TO/THE/NisslST/DIRECTORY';
addpath(genpath(nissl_st_dir));
output_dir = fullfile(nissl_st_dir,'output');
tractography_dir = fullfile(output_dir,'tractography');
if ~exist(tractography_dir,'dir')
    mkdir(tractography_dir)
end


%% Load the structure tensor statistics file
ds_factor = 1;% Downsample factor for downsampling the image
pix_microns = 0.97;
side_microns = 200; % Tile side in microns
rho = 15; % measured in pixels (~15 microns)
sigma = 3; % Gaussian kernel for smoothing the image before calculating Nissl-ST. This parameter is typically 0, but for this Allen Brain Atlas dataset, it seems that sigma=3 yields better results.
nPeaks = 2; % Number of peak orientations to extract. In fact, we typically just take the first peak.
sample_near_cells = true; % Only extract orientations as close as 'rho' to cells (to avoid measuring empty spaces)

st_stats_file = fullfile(output_dir,sprintf('st_stats_side_%gmu_rho%g_sigma%g_npeaks%g_samplenearcells_%g_slicethresh_ds%g.mat',side_microns,rho,sigma,nPeaks,sample_near_cells,ds_factor));
if ~exist(st_stats_file,'file')
    error('File does not exist. Run the NisslST_run_exampl.m first')
end
load(st_stats_file);

%% Get initial brain mask based on the Otsu threshold per sub-image
% Note: The functions used for creating the masks typically work best when
% using entire slices.
mask = getMaskFromThreshMap(thresh_map);

% Get GM mask based on the distribution of cluster cell sizes in pixels (look at 95th prctile)
cluster_pixel_size_map = gm_wm_info.cluster_pixel_size(:,:,4);
[gmMask, wmMask] = get_gm_wm_masks(cluster_pixel_size_map,mask);

% Get a grayscale image with GM/WM contrast
grayim = squeeze(cluster_pixel_size_map(:,:,1));
grayim(abs(grayim)>400) = 0;
grayim(~mask) = 0;

%% Run tractography on the peak orientations
% This requires the vistasoft toolbox for writing NIfTI files
% (https://github.com/vistalab/vistasoft). There are probably alternatives
% that could be just as good.

% Export a NIfTI file of the peak orientations
% This can be used for tractography using MRTrix's implementation of FACT.
fact_input_file = fullfile(tractography_dir,'FACT_input.nii.gz');
write_FACT_input(theta_peaks_map,wmMask,fact_input_file);

% Write the wmMask to file
mask_file = fullfile(tractography_dir,'tractography_mask.nii.gz');
xform = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]; % Because these are coronal sections, treat the 2nd dimension as z axis
dtiWriteNiftiWrapper(wmMask,xform,mask_file);

% Write a seed mask to file
% NOTE: A seed mask is typically drawn manually, but here we specify it for
% this example
seed_mask = zeros(size(wmMask));
% seed_mask(15,33:42) = 1;
seed_mask(31,59:69) = 1;
xform = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]; % Because these are coronal sections, treat the 2nd dimension as z axis
seed_mask_file = fullfile(tractography_dir,'seed_mask.nii.gz');
dtiWriteNiftiWrapper(seed_mask,xform,seed_mask_file);


% Write the gmMask to file
gmMask_file = fullfile(tractography_dir,'gray_matter_mask.nii.gz');
dtiWriteNiftiWrapper(double(imdilate(gmMask,strel('disk',0))),xform,gmMask_file);

% Write a system command for running tckgen
% ****
seed_name = 'seed_mask'; %'tractography_mask'; % If you haven't created a seed mask, use the white-matter mask as default
seed_file = fullfile(tractography_dir,[seed_name, '.nii.gz']);
cfg = [];
cfg.maxlength = '5000';
cfg.step = '3';
cfg.angle = '60';
cfg.seed_grid_per_voxel = [seed_file ' 5'];
cfg.select = '10';
cfg.mask = mask_file;
% cfg.include =
% {fullfile(tractography_dir,sliceName,'waypoint_midT_supinf.nii.gz')}; % Only retain streamlines that pass through this ROI
% cfg.include = {gmMask_file}; % Only retain streamlines that reach the gray matter mask. In this case of a small image that does not include the entire slice, it makes no sense to require that
tck_out = fullfile(tractography_dir,[seed_name '_step' cfg.step '_angle' cfg.angle '.tck']);
cmd = st_write_tck_cmd(fact_input_file,tck_out,cfg);
% If working in Windoes:
cmd = strrep(cmd,'\','\\');
cmd = ['C:\\msys64\\mingw64\\bin\\' cmd]; 
system(cmd) 
% NOTE: If this doesn't work, you might need to install MRTrix 3.0 through Msys64 and paste the command over there
% NOTE: If you the tckgen command is not recognized, try using the full path of the command, for example: cmd = ['C:\\msys64\\mingw64\\bin\\' cmd];


%% Load the tractography results
fg = fgRead(tck_out);
%'C:\Users\roey\Documents\MezerLab\Glia_orientation_paper\Code_for_github\NisslST-main\output\tractography\seed_roi_step3_angle60.tck');

% Flip left-right to get correct orientation:
fg.fibers = cellfun(@(x) bsxfun(@times,x,[-1,1,1]'), fg.fibers,'UniformOutput',false);

% fg.fibers = fg.fibers(1:100:end); % Choose subset to plot
%% Get RGB value for each node, according to its local orientation
rgb = [];
for fI = 1:length(fg.fibers)
   fib = fg.fibers{fI};
%    fib(1,:) = -fib(1,:); % Flip the x axis to get the correct colors
   theta_vec = []; 
   for nI = 1:length(fib)-1
        tmp = fib(:,nI+1)-fib(:,nI);
        tmp(2) = [];
         theta = atan2(tmp(2),tmp(1))*180/pi;
         if theta<0
            theta = theta+180; 
         end
         vals = [theta; 0; 180];
         rgbTmp = vals2colormap(vals, hsv, []);
         rgbTmp(end-1:end,:) = [];
    rgb{fI}(nI,:) = rgbTmp;
   end
    rgb{fI}(end+1,:) = rgb{fI}(end,:);
end

%% Plot the streamlines
fg.fibers = cellfun(@(x) bsxfun(@times,x,[-1,1,-1]'), fg.fibers,'UniformOutput',false);
hlight = AFQ_RenderFibers(fg,'color',rgb,'tubes', 1, 'jittershading',0,'radius',[0.3,1],'camera','coronal', 'newfig',1);
set(gcf,'position',get(0,'screensize'));
grid off; axis off
set(gca,'color','k')
set(gcf,'color','k')
