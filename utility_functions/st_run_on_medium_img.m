% THIS CODE RUNS THROUGH EACH SUBDIRECTORY (A LARGE SUBIMAGE) OF SOME SLICE. IT THEN EXTRACTS
% GLIA ORIENTATION USING THE STRUCTURE TENSOR FOR EACH SUB-IMAGE AND SAVES INFO FOR ALL SUB-IMAGES
% TO A FILE. LATER ON, THE CODE THAT PLOTS THE ORIENTATIONS MAPS LOADS THIS FILE AND AN ROI MASK TO PLOT % % ONLY RELEVANT SUB-IMAGES.
% clear; clc
function st_run_on_medium_img(imFile,pix_microns,side_microns,rho,sigma,nPeaks,sample_near_cells,ds_factor,st_stats_file)
% dI is which diretory to run (each directory is a slice)
% imI_toRun is an index (or indices) between 1 an 50

%% Convert microns to pixels
side_pix = round(side_microns/pix_microns); % The side of each sub-image in pixels

%% Find directy do work on
% If this is a downsampled image, modify the parameters
if ds_factor>1
    side_pix = round(side_microns/(pix_microns/ds_factor)); % The side of each sub-image in pixels
    rho = round(rho*ds_factor);
end

% Extract structure tensor statistics
[theta_mean_map, aniso_map, coherence_map, theta_peaks_map, peaks_height_map, theta_von_mises_map, von_mises_component_fraction_map, thresh_map, mean_val_map] = st_get_img_structure_tensor_statistics(imFile,side_pix,rho,sigma,nPeaks,sample_near_cells,ds_factor);
% Extract informationn to get a grey-matter mask
gm_wm_info = st_get_img_gm_wm_info(imFile,side_pix);

%% Save
save(st_stats_file,'theta_mean_map','aniso_map','coherence_map','theta_peaks_map','peaks_height_map','thresh_map','theta_von_mises_map','von_mises_component_fraction_map','gm_wm_info','mean_val_map');
    
end
