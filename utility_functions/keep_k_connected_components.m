function mask = keep_k_connected_components(mask,k)
CC = bwconncomp(mask);
k = min(k,numel(CC.PixelIdxList));% I case there are less than k connected components
len = cellfun(@length, CC.PixelIdxList);
[~,idx] = sort(len,'descend');
mask = zeros(size(mask));
for cI = 1:k
    mask(CC.PixelIdxList{idx(cI)}) = 1;
end
mask = double(logical(mask));
end