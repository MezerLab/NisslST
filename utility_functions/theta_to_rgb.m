function rgb = theta_to_rgb(theta_map,mask)
% Input is nxmx1 or in general nxmxp matrix of angles in the range [0,180]
% Output is an rgb image, masked to the brain slice
colormapdim = 2; % 2D colormap (hsv, as uauall), or 3D map (blue-red for coronal slices)
rgb = [];
for vI = 1:size(theta_map,3)
    vals = theta_map(:,:,vI);
    vals = vals(:);
    vals = [vals; 0; 180];
    if colormapdim==3
        xx = cosd(vals);
        zz = sind(vals);
        xx = abs(xx);
        zz = abs(zz);
        rgbTmp = [1 0 0].*xx + [0 0 1].*zz;
        rgbTmp(end-1:end,:) = [];
    elseif colormapdim==2
        rgbTmp = vals2colormap(vals, hsv, []);
        rgbTmp(end-1:end,:) = [];
    end
    
    rgbTmp = reshape(rgbTmp,[size(theta_map(:,:,vI),1),size(theta_map(:,:,vI),2),3]);
    [ii,jj] = find(~mask);
    for idx = 1:length(ii)
        rgbTmp(ii(idx),jj(idx),:) = [0 0 0];
    end
    
    % Remove exactly 90 degrees
    %     [ii,jj] = find(theta_map(:,:,vI)==90);
    %     for idx = 1:length(ii)
    %         rgbTmp(ii(idx),jj(idx),:) = [0 0 0];
    %     end
    
    rgb{vI} = rgbTmp;
end
end