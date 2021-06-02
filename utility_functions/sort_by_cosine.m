function theta_map = sort_by_cosine(theta_map)
% Sort an array of mXnX2 of angles is [0,180] by their cosine, to separate
% high and low angles.
theta_map_cos = cosd(theta_map);
tmp = theta_map;
for d1 = 1:size(theta_map,1)
    for d2 = 1:size(theta_map,2)
        if theta_map_cos(d1,d2,1) < theta_map_cos(d1,d2,2)
            tmp(d1,d2,1) = theta_map(d1,d2,2);
            tmp(d1,d2,2) = theta_map(d1,d2,1);
        end
    end
end
theta_map = tmp;
end

