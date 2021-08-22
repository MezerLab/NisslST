function im = replace_holes_with_nearest_neighbor(im,mask)
[x_all,y_all] = find((isnan(im) & mask) | im==0 & mask); % pixels to replace
for pI = 1:length(x_all)
   x = x_all(pI);
   y = y_all(pI);
    vals = [im(x-1,y-1),im(x-1,y),im(x-1,y+1), im(x,y-1),im(x,y+1),im(x+1,y-1),im(x+1,y),im(x+1,y+1)];
   vals(isnan(vals)) = [];
   im(x,y) = mode(vals);
end
end