function st_run_on_medium_img(imFile,pix_microns,side_microns,rho,sigma,nPeaks,sample_near_cells,ds_factor,st_stats_file)
% ST_RUN_ON_MEDIUM_IMG performs Nissl Structure-Tensor (Nissl-ST) analysis
% on a an input image. It runs separately on evey image tile (typically
% 50-200^2 microns^2) in the input image.
%
% INPUT
% =====
% imFile: The input image. Can be a file or an image matrix (e.g., the
%         output of MATLAB's imread)
% pix_microns: The raw-data resolution. Pixel size in Atlas of the Human
%              Brain is 0.645 microns. In Allen Brain Atlas, it is typically 0.97
%              microns.
% side_microns: The side of the image tiles to you (typically 50-200
%               microns)
% rho: The standard deviation of the Gaussian kernel used for spatial
%      regularization in the pixel-wise structure tensor calculation, in
%      pixels. Typically 15 microns, so for data with 0.645 micron
%      resolution, use rho = 23 means 23*0.645 =~ 15 mu.
% sigma: The standard deviation of the Gaussian kernel used for smoothing
%        the image prior to calculating the structure tensor. Typically 0.
% nPeaks: The number of peak orientations to extract. We typically enter 2,
%         but actually only look at the first peak. Exploring additional
%         peaks requires further investigation.
% sample_near_cells: Restrict the analysis to pixels within a distance of
%                    rho from any stained cell. This allows ignoring pixels
%                    with only background tissue, which might bias the
%                    estimated peak orientation in the tile.
% ds_factor: Downsampling factor. Downsample the input image by this factor
%            to explore the effect of input resolution on Nissl-ST.
%            Default: 1 (no downsampling)
% st_stats_file: Full path of the desired output .mat file.
%
% OUTPUT 
% ====== 
% The function saves the analysis results in st_stats_file, which includes: 
% theta_peaks_map: An N x M x nPeaks map with peak orientations for each of
%                  the NxM image tiles (in the range of [0,180] degrees.
% coherence_map: Coherence of the pixelwise orientations in each tile.
%                Coherence is the norm of the vector sum of all pixelwise
%                orientations. It ranges between 0 (incoherent) and 1
%                (coherent orientations).
% gm_wm_info: A structure with maps that are used for calculating the
%             white-matter and gray-matter masks. Sepcifically, we binarize
%             each image tile, run connected components analysis, and
%             extract a measure of the size of the identified objects. In
%             the gray-matter, where a larger area of the cells are
%             stained, and nearyby cells tends to form larger connected
%             objects in the image, this value tends to be high. Therefore,
%             plotting the map of object sizes results in an image with a
%             high contrast between gray- and white-matter, which
%             fascilitates the creation of a white-matter mask.
% Additional outputs are rarely used, but might be helpful in further developing the analysis:
% peaks_height_map: The normalized height of each identified peak in each
%                   tile. This allows ordering the peaks by prodeominance.
%                   Notice that other measures might be better for this
%                   purpose. For example, the peak height multiplied by the
%                   peak width, or estimated area under the curve for each
%                   peak, if modelled as a  mixture model.
% theta_mean_map: The orientation of the vector-sum of all pixelwise
%                 orientations in each image tile.
% aniso_map: The average local anisotropy in each image tile. This is
%            defined as the ratio of the eigenvalues of the structure
%            tensor.
% thresh_map: The Otsu threshold that optimizes the binarization of each
%             image tile. When plotted, gives some contrast between white-
%             and gray-matter.
% theta_von_mises_map: Similar to theta_peaks_map, but with peak
%                      orientations identifies as a mixture of von Mises
%                      distributions. This is still under development.
% von_mises_component_fraction_map: The estimated component fraction of
%                                   each component in the mixture model.
% mean_val_map: The mean gray-scale value in each image tile. When plotted,
%               gives some contrast between white- and gray-matter.
% See st_get_img_structure_tensor_statistics.

%% Convert microns to pixels
side_pix = round(side_microns/pix_microns); % The side of each sub-image in pixels

%% Perform Nissl-ST
% Extract structure tensor statistics
[theta_mean_map, aniso_map, coherence_map, theta_peaks_map, peaks_height_map, theta_von_mises_map, von_mises_component_fraction_map, thresh_map, mean_val_map] = st_get_img_structure_tensor_statistics(imFile,side_pix,rho,sigma,nPeaks,sample_near_cells,ds_factor);
% Extract informationn to get a grey-matter mask
gm_wm_info = st_get_img_gm_wm_info(imFile,side_pix);

%% Save
save(st_stats_file,'theta_mean_map','aniso_map','coherence_map','theta_peaks_map','peaks_height_map','thresh_map','theta_von_mises_map','von_mises_component_fraction_map','gm_wm_info','mean_val_map');
    
end
