function img = reorient_img(img)
tmp = [];
for dI = 1:size(img,3)
    tmp(:,:,dI) = imrotate(squeeze(img(:,:,dI)),-90);
    tmp(:,:,dI) = flipdim(squeeze(tmp(:,:,dI)),2);
end
img = tmp;
end
