
%% Another version of this function, that takes as input a large image instead of a directory with subimages
function [theta_mean_map, aniso_map, coherence_map, theta_peaks_map, peaks_height_map, theta_von_mises_map, component_fraction_map, thresh_map, mean_val_map] = st_get_img_structure_tensor_statistics(imFile,side_pix,rho,sigma,nPeaks,sample_near_cells,ds_factor,varargin)
% st_get_img_structure_tensor_statistics performs Nissl-ST and returns
% multiple outputs. See st_run_on_mediam_img for details on all outputs.
%
% If varargin>0, then this function also plots test plots to investigate
% the output, including a polar histogram of the pixelwise orientations.

%%
% See if this is a special call to the function only to make test plots 
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
im = im_gray;

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

if ds_factor~=1 % If required to downsample the image tile fix the rho parameter
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
            try
                [~,fig_name] = fileparts(imFile);
            catch
                fig_name = 'test';
            end
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
        
        % Plot circular histogram (this requires MATLAB version 2016b and
        % above)
        if call_plot_version 
            if exist('polarhistogram','builtin')
            figure('color','w')
            theta_plot = theta_vec(:);
            theta_plot = [theta_plot; theta_plot(theta_plot<0)+180; theta_plot(theta_plot>0)-180];
            theta_plot = theta_plot.*pi/180;
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
            set(gca,'rlim',[0,prod(size(sib_img))/8]); % 310^2 is number of pixels in the image. Change this if using a different image tile size, or a different dataset.
            set(gca,'color','k')
            h.EdgeColor = 'w';
            h.FaceColor = 'w';
        
            % Add the peak orientations
            for mI = 1
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
        else
            warning('The function ''polarhistogram '' is required for plotting the polar histogram, but this version of MATLAB does not have it. Skipping plot')
            end
        end
    end
end

end

%% Functions
function [theta_peaks,pks_height] = find_peak_orientations(theta_vec,nPeaks)
% Input is in [0 180], so we expand it below. 
    % Convert to radians
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