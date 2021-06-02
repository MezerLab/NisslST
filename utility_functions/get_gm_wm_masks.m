function [gmMask, wmMask] = get_gm_wm_masks(cluster_pixel_size_map,mask)
cluster_pixel_size_map(cluster_pixel_size_map>1000) = 0; % Remove outliers
gmMask = cluster_pixel_size_map(:,:,end).*mask;
tmp = gmMask(gmMask>0);
% Find threshol to separte GM from WM
bins_centers = [1:254]+0.5;
counts = hist(tmp(:),bins_centers);
[gm_thresh,~] = otsuthresh(counts);
gm_thresh = gm_thresh.*255;
% Create a WM mask
wmMask = zeros(size(gmMask));
wmMask(gmMask<gm_thresh & gmMask>0) = 1;
% Create a GM mask
tmp = zeros(size(gmMask));
tmp(gmMask>gm_thresh) = 1;
gmMask = tmp;

% Clean the edges (maybe we can do with a little less cleaning)
% wmMask= imerode(wmMask,strel('disk',1));
% wmMask= imdilate(wmMask,strel('disk',1));
wmMask= imfill(wmMask,'holes');
k = 10;
wmMask = keep_k_connected_components(wmMask,k);
%figure
%imagesc(wmMask)

% Finally, remove WM from the GM
gmMask(find(wmMask)) = 0;
%%
end



