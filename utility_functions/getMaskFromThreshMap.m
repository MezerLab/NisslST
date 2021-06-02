function mask = getMaskFromThreshMap(thresh_map)
% This creates a brain mask
mask = thresh_map;
mask = imbinarize(mask,'global');

% Erode and dilate the image to eliminate small noisy regions
erode_dilate_num = 3;
se = strel('disk',1);
for ii = 1:erode_dilate_num
    mask = imerode(mask,se);
end
for ii = 1:erode_dilate_num
    mask = imdilate(mask,se);
end

% Remove small disconnected components
cc = bwconncomp(mask);
len = cellfun(@length,cc.PixelIdxList);
small_clusters = find(len<300); % This threshold relly depends on the size of the image you're working with. Look at the distribution of sizes of the connected components you get to figure out what size should be considered as noise
for cI = 1:length(small_clusters)
    mask(cc.PixelIdxList{small_clusters(cI)}) = 0;
end

% Clean the edges (maybe we can do with a little less cleaning)
mask = imerode(mask,strel('disk',2));

% Fill holes
mask = imfill(mask,'holes');
end
