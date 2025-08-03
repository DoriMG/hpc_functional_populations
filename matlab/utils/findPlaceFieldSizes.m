function [sizePf] = findPlaceFieldSizes(pcLoc, bin_size)
% Adapted from Bourboulou et al. (2019)
% (https://doi.org/10.7554/eLife.44487)
sizePf = nan(size(pcLoc,1),1);
for i = 1:size(pcLoc,1)
    data_smooth = smooth(pcLoc(i,:),0.13,'loess')'; % A wee smooth
    if sum(isnan(data_smooth)) < .5*length(data_smooth)
        % Find the half max value.
        halfMax = (min(data_smooth) + max(data_smooth)) / 2;
        [~,centre] = nanmax(data_smooth);
        
        %find where the data smooth line (with more points) first crosses the half
        %max line and last crosses the line
        [~, idxs] = find(data_smooth<=halfMax);
        
        % If the place field touches teh edge of the env, set the edge of 
        %env as the outside edge of teh place field
        outside_idxs = idxs-centre;
        if sum(outside_idxs<0) > 0
            index1 = centre + max(outside_idxs(outside_idxs<0));
        else
            index1 = 1;
        end
        
        if sum(outside_idxs>0) >0
            index2 = centre + min(outside_idxs(outside_idxs>0));
        else
            index2 = length(data_smooth)+1;
        end
        
        % Calculate the width
        width = (index2-index1)*bin_size;
        sizePf(i) = width;
    end
end
end

