
%% Another version of this function, that takes as input a large image instead of a directory with subimages
function [gm_wm_info] = st_get_img_gm_wm_info(imFile,side_pix)
% imFile is the full path of an image file to be cropped into sub-images,rh
% side_pix is the length of the sub-image to be cropped in pixels (should
% be 155 for the AHB dataset, to get 50 mu sub-images)

if isstr(imFile)
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

% Binarize the image to get cells masks
% % % im = imbinarize(im_gray,'global'); &&& This is what I had before, but
% gave spatially varying results
thresh = 0.7; 
im = im2bw(im_gray,gray(255),thresh); % &&& Using same thresh for all images 
im = 1-im; % Cells are now bright

%### Roey End of testing contrast enhancement
side_pix = side_pix-1; % The imcrop function adds 1 by default

iSize = floor(size(im,2)/side_pix); % x axis is the 2nd dimension (columns)
jSize = floor(size(im,1)/side_pix);

lenMat = [];

for ii = 1:iSize
    for jj = 1:jSize
        upper_left_x = 1+side_pix*(ii-1);
        upper_left_y = 1+side_pix*(jj-1);
        rect = [upper_left_x, upper_left_y, side_pix, side_pix];
        imTmp = imcrop(im,rect);
        
        % Connected components
        CC = bwconncomp(imTmp);
        
        len = cellfun(@length, CC.PixelIdxList); % Size of each patch in pixels
        lenMat(ii,jj,1) = mean(len);
        lenMat(ii,jj,2:4) = prctile(len,[80,90,95]);
    end
end
%     if ~mod(fI,100)
%         orientation_and_aniso_file = fullfile(imDir,'ahb_r3-0303_orientation_and_anisotropy.mat');
%         save(orientation_and_aniso_file,'rgb_orientation_map','aniso_map');
%     end
gm_wm_info.cluster_pixel_size = lenMat;

end
