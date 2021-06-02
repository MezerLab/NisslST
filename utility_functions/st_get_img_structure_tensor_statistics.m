
%% Another version of this function, that takes as input a large image instead of a directory with subimages
function [theta_mean_map, aniso_map, coherence_map, theta_peaks_map, peaks_height_map, theta_von_mises_map, component_fraction_map, thresh_map, mean_val_map] = st_get_img_structure_tensor_statistics(imFile,side_pix,rho,sigma,nPeaks,sample_near_cells,ds_factor,varargin)
% imFile is the full path of an image file to be cropped into sub-images,rh
% side_pix is the length of the sub-image to be cropped in pixels (should
% be 155 for the AHB dataset, to get 50 mu sub-images)
% rho is the radius(?) of the environment which is used for calculating the
% local structure tensors. rho = 23 means 23*0.645 =~ 15 mu

% See if this is a special call to the funciton only to make plots for the methods section
call_plot_version = false;
if numel(varargin)>0
   call_plot_version = true; 
end

if ischar(imFile)
    im = imread(imFile);
else
    im = imFile; %If an image was given as input
end

if length(size(im))==3 % im is rgb
    im_gray = rgb2gray(im);
else
    im_gray = im;    
end
% im_gray  = imadjust(im_gray); % imadjust
% im_gray  = histeq(im_gray); % histogram equalization
im = im_gray;



%### Roey End of testing contrast enhancement
side_pix = side_pix-1; % The imcrop function adds 1 by default

iSize = floor(size(im,2)/side_pix); % x axis is the 2nd dimension (columns)
jSize = floor(size(im,1)/side_pix);

% Initialize
mean_val_map = nan(iSize,jSize);
thresh_map = nan(iSize,jSize);
theta_mean_map = nan(iSize,jSize);
aniso_map= nan(iSize,jSize);
coherence_map= nan(iSize,jSize);
theta_peaks_map = nan(iSize,jSize,nPeaks);
peaks_height_map = nan(iSize,jSize,nPeaks);
theta_von_mises_map = nan(iSize,jSize,nPeaks);
component_fraction_map = nan(iSize,jSize,nPeaks);

if ds_factor~=1 % If required to downsample the sub-image, fix the rho parameter
    rho = round(rho*ds_factor);
end

for ii = 1:iSize
    for jj = 1:jSize
        upper_left_x = 1+side_pix*(ii-1);
        upper_left_y = 1+side_pix*(jj-1);
        rect = [upper_left_x, upper_left_y, side_pix, side_pix];
        sub_img = imcrop(im,rect);
        % Skip empty images
        if length(unique(sub_img(:)))==1
            continue
        end

        % Denoise for abatlas &&&
%         sub_img = BM3D(sub_img, 5);
%         sub_img = sub_img./max(sub_img(:));

        % Downsample the sub-image, if required
        if ds_factor~=1
            sub_img = imresize(sub_img,ds_factor);
        end
        
        
        % Extract structure tensor statistics
        if ~call_plot_version
            % For regular analysis
            [theta_mean, aniso_mean, coherence, theta_vec, thresh, mean_val] = st_get_subimg_structure_tensor_statistics(sub_img,rho,sigma,sample_near_cells);
        else
            % For getting figures:
            % Name of this sub-figure
            [~,fig_name] = fileparts(imFile);
            fig_name = sprintf([fig_name,'subidx_ii%g_jj%g_rho%g_sigma_%g_sample_near_cells%g'],ii,jj,rho,sigma,sample_near_cells);
            [theta_mean, aniso_mean, coherence, theta_vec, thresh, mean_val] = st_get_subimg_structure_tensor_statistics_figures_version(sub_img,rho,sigma,sample_near_cells,fig_name);
        end
        
        if isempty(theta_vec) || all(isnan(theta_vec))
            continue
        end
        % Extract the peak orientation (from the possibly masked eigvec_weighted) 
        [theta_peaks,pks_height] = find_peak_orientations(theta_vec,nPeaks);
        % Fit von Mises distribution to the eigvecs
        clear vonMises
        try
            theta_vec(theta_vec<90)=theta_vec(theta_vec<90)-180; % Get it back to [-pi,pi] instead of [0,pi]
            vonMises = fitmvmdist(theta_vec*pi/180,nPeaks,'MaxIter',250); % Set maximum number of EM iterations to 250
        catch
            vonMises.componentProportion = zeros(1,nPeaks);
            vonMises.mu = zeros(1,nPeaks);
        end
        [componentProportion,I] = sort(vonMises.componentProportion,'descend');
        vonMises_peaks = vonMises.mu(I);
        
        % Insert all statistics into maps
        theta_mean_map(ii,jj) = theta_mean;
        aniso_map(ii,jj) = aniso_mean;
        coherence_map(ii,jj) = coherence;
        theta_peaks_map(ii,jj,:) = theta_peaks; % Ordered by peak height
        peaks_height_map(ii,jj,:) = pks_height;
        theta_von_mises_map(ii,jj,:) = vonMises_peaks*180/pi; % Ordered by proportion (may be different than height, e.g. for a very broad and low-peak distribution)
        component_fraction_map(ii,jj,:) = componentProportion;
        thresh_map(ii,jj) = thresh;
        mean_val_map(ii,jj) = mean_val;
        
        % Plot circular histogram
        plot_circ_hist = false;
        if plot_circ_hist | call_plot_version
            figDir = '/ems/elsc-labs/mezer-a/Mezer-Lab/analysis/Histology/Figures_method/output_figures';
            figDir = '/ems/elsc-labs/mezer-a/Mezer-Lab/analysis/Histology/Figures_method/parameter_figure/';
%             figDir = '/ems/elsc-labs/mezer-a/Mezer-Lab/analysis/Histology/Figures_method/parameter_figure_bw/';

            rho_sigma_str = sprintf('rho%g_sigma%g',rho,sigma);
            %%
            figure('color','w')
            theta_plot = theta_vec(:);
            theta_plot = [theta_plot; theta_plot(theta_plot<0)+180; theta_plot(theta_plot>0)-180];
            theta_plot = theta_plot.*pi/180;
%             circ_plot(theta_plot,'hist',[],60,true,true,'linewidth',2,'color','r')
            bins_num = 60;
            bins_edges = [0:3:360]./180*pi;
            h = polarhistogram(theta_plot,'BinEdges',bins_edges,'Normalization','count','linew',3);
            h.DisplayStyle = 'bar'; % Or 'stairs' for unfilled histogram
            h.EdgeColor = [0,0,0];
            h.FaceColor = [0,0,0];
            h.FaceAlpha = 1;
            set(gca,'RTick',[])
            set(gca,'ThetaTickLabel',[])
            set(gca,'ThetaMinorTick','on')
            set(gca,'linew',2)
            hold all
            mu_deg = theta_peaks;
            rl = rlim;
            set(gca,'rlim',[0,(310^2)/8]); % 310^2 is number of pixels in the image &&& CHange for different datasets
            set(gca,'color','k')
            h.EdgeColor = 'w';
            h.FaceColor = 'w';
            export_fig(gcf,fullfile(figDir,[fig_name,'_',rho_sigma_str,'_histogram_sample_near_cells.png']),'-dpng','-r300');
        
%             h.DisplayStyle = 'stairs'; % Or 'stairs' for unfilled histogram
%             cmap = lines(8);
%             h.EdgeColor = cmap(4,:); % 3 for grayscale, 4 for binary 
%             export_fig(gcf,fullfile(figDir,[fig_name,'_',rho_sigma_str,'_histogram_stairs_count_sample_near_cells.png']),'-dpng','-r300');
            close(gcf) %&&&DELETE ME if contiuing to add peaks later
            continue
        
            % Add the peak orientations
            for mI = 1:2
                if mI == 1;frmt = '-r';else;frmt = '--r'; end
                hplt = polarplot([mu_deg(mI),mu_deg(mI)-180]*pi/180,[1,1].*rl(2),frmt,'linew',3)
                % Get orientation RGB color
                vals = mu_deg(mI);
                vals = [vals, 0, 180];
                rgbTmp = vals2colormap(vals, hsv, []);
                rgbTmp(end-1:end,:) = [];
                hplt.Color = rgbTmp;
                %plot([0, cosd(mu_deg(mI))],[0, sind(mu_deg(mI))],frmt,'linew',2)
                set(gca,'color','w')
                h.EdgeColor = 'k';
                h.FaceColor = 'k';
                export_fig(gcf,fullfile(figDir,[fig_name,'_',rho_sigma_str,'_histogram_peaks_',num2str(mI),'_sample_near_cells.png']),'-dpng','-r300');
                set(gca,'color','k')
                h.EdgeColor = 'w';
                h.FaceColor = 'w';
                export_fig(gcf,fullfile(figDir,[fig_name,'_',rho_sigma_str,'_histogram_peaks_',num2str(mI),'_black_sample_near_cells.png']),'-dpng','-r300');
            end
            %%
            close gcf
%             mu_deg = fittedVmm.mu*180/pi;
%             [~,idx]= sort(fittedVmm.componentProportion,'descend');
%             for mI = idx'
%                 if mI == 1;frmt = '-k';else;frmt = '--k'; end
%                 plot([0, cosd(mu_deg(mI))],[0, sind(mu_deg(mI))],frmt,'linew',2)
%             end
        end
    end
end

end

%% Functions
function [theta_peaks,pks_height] = find_peak_orientations(theta_vec,nPeaks)
% Input is in [0 180], so we expand it below. % Convert to radians
    theta_vec = theta_vec/180*pi;
    % Apply ksdensity to smooth the histogram and find the optimal
    % smoothing kernel bandwidth
    [~,~,bw] = ksdensity(theta_vec,0:pi/180:pi); % bw is the smoothing kernel bandwidth
    % Initialize
    theta_peaks = nan(1,nPeaks);
    pks_height = nan(1,nPeaks);
    % Use ksdensity to smooth histogram and find peaks
    theta_tmp = [theta_vec-pi; theta_vec; theta_vec+pi]; % Repeat the data, so we don't interpret peaks near 0 as two peaks
    xi = -pi : pi/180 : 2*pi;
    [f] = ksdensity(theta_tmp,xi,'BandWidth',bw);
        
    % Normalize to maximal peak of 1
    f = f./max(f);
    % Find peaks
    [pks,locs] = findpeaks(f);
    if isempty(pks)
        return
    end
    % Remove peaks outside of [0,pi]
    idx = find(xi(locs)>0 & xi(locs)<pi);
    locs = locs(idx);
    pks = f(locs);
    % Sort peaks in descending order (consider ordering according to width
    % multiplied by peak, to estimate the integral under the peak). See
    % findpeaks for the width.
    [pks,I] = sort(pks,'descend');
    locs = locs(I);
    for pI = 1:length(locs)
        theta_peaks(pI) = xi(locs(pI));
        pks_height(pI) = pks(pI);
    end
    pks_height = pks_height (1:nPeaks);
    theta_peaks = theta_peaks(1:nPeaks);
    theta_peaks = theta_peaks*180/pi;
    %%
    test_plots = false;
    if test_plots
        figure('color','w')
        histogram(theta_tmp);
        [f] = ksdensity(theta_tmp,xi,'BandWidth',bw);
        hold all
        scale = get(gca,'ylim');
        plot(xi,f*scale(2)/2)
        for pI = 1:length(theta_peaks)
           if isnan(theta_peaks(pI))
               continue
           end
           plot(theta_peaks(pI)*[1 1],ylim);
        end
    end
    %%
end