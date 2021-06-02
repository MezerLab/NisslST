function imshow_rgb_on_gray(rgbim,grayim,mask,varargin)
% Plot rgb image on top of gray image, where mask is where rgb image should
% be shown
imshow(rgbim)
if numel(varargin)>1
    clims = varargin{2};
    caxis(clims);
end
hold on
% grayim = floor(grayim./max(grayim(:))*255);
% grayim = cat(3, grayim, grayim, grayim);
h = imshow(grayim);
vals = grayim(mask);
if numel(varargin)>0
    clims = varargin{1};
else
    clims = prctile(vals,[5,99.9]);
end
caxis(clims);
hold off
pause(0.25);
set(h,'AlphaData',1-mask);
end