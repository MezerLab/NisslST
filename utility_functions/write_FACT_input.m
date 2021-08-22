function write_FACT_input(theta_peaks_map,wmMask,fact_input_file)
% WRITE_FACT_INPUT takes as input a map of orientation peaks and a white-matter mask, 
% and writes a volume nifti image to be used as input to FACT tractography.
% see: https://mrtrix.readthedocs.io/en/latest/reference/commands/tckgen.html

%% The old version
% % % % Fill holes with nearest neighbor interpolation
% % % theta_peaks_map(:,:,1) = replace_holes_with_nearest_neighbor(theta_peaks_map(:,:,1),wmMask);
% % % xform = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]; % Becasue these are coronal sections, treat the 2nd dimension as z-axis
% % % 
% % % sz = size(theta_peaks_map);
% % % tmp = zeros(sz(1),sz(2),1,3);
% % % for ii  = 1:sz(1)
% % %     for jj = 1:sz(2)
% % %         ang = theta_peaks_map(ii,jj,1);
% % %         if isnan(ang)
% % %             tmp(ii,jj,:,1) = 0;
% % %             tmp(ii,jj,:,2) = 0;
% % %             tmp(ii,jj,:,3) = 0;
% % %             continue
% % %         end
% % % mat = rotz(-ang)*[1 0 0]';
% % % tmp(ii,jj,:,1) = mat(1);
% % % tmp(ii,jj,:,2) = mat(3); % 2&3 dimensions are exchanged, because this is a coronal plane ***
% % % tmp(ii,jj,:,3) = mat(2);
% % %     end
% % % end
% % % 
% % % tmp = bsxfun(@times,tmp,wmMask);
% % % % tmp = tmp/1000;
% % % dtiWriteNiftiWrapper(double(tmp),xform,fact_input_file);

%% The new version I'm trying to get tractography to work properly
% Fill holes with nearest neighbor interpolation
theta_peaks_map(:,:,1) = replace_holes_with_nearest_neighbor(theta_peaks_map(:,:,1),wmMask);
xform = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]; % Becasue these are coronal sections, treat the 2nd dimension as z-axis
% xform = eye(4);

sz = size(theta_peaks_map);
tmp = zeros(sz(1),sz(2),1,3);
for ii  = 1:sz(1)
    for jj = 1:sz(2)
        ang = theta_peaks_map(ii,jj,1);
        if isnan(ang)
            tmp(ii,jj,:,1) = 0;
            tmp(ii,jj,:,2) = 0;
            tmp(ii,jj,:,3) = 0;
            continue
        end
mat = roty(ang)*[1 0 0]'; %% Use roty because this is a coronal slice
tmp(ii,jj,:,1) = mat(1);
tmp(ii,jj,:,2) = mat(2); % 2&3 dimensions are exchanged, because this is a coronal plane ***
tmp(ii,jj,:,3) = mat(3);
    end
end

tmp = bsxfun(@times,tmp,wmMask);
% tmp = tmp/1000;
dtiWriteNiftiWrapper(double(tmp),xform,fact_input_file);
end