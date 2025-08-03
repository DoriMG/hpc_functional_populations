function [inOutRate] = findOutInFieldRate(pcLoc)
% Adapted from Bourboulou et al. (2019)
% (https://doi.org/10.7554/eLife.44487)
inOutRate = nan(size(pcLoc,1),1);
for i = 1:size(pcLoc,1)
    data_smooth = smooth(pcLoc(i,:),0.13,'loess')';
    
    if sum(isnan(data_smooth)) < .5*length(data_smooth)
        % Find the half max value.
        halfMax = (min(data_smooth) + max(data_smooth)) / 2;
        [~,centre] = nanmax(data_smooth);
        %find where the data smooth line (with more points) first crosses the half
        %max line and last crosses the line
        [~, idxs] = find(data_smooth<=halfMax);
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
        
        % find what bins are in the field
        inField = zeros(length(data_smooth),1);
        inField(index1:(index2-1)) = 1;
        
        % calculate the normalized in/out ratio per cell
        data_norm = pcLoc(i,:) - min(pcLoc(i,:));
        data_norm = data_norm./max(data_norm);
        outAct = nanmean(data_norm(~logical(inField)));
        inAct = nanmean(data_norm(logical(inField)));
        inOutRate(i) = outAct/inAct;
    end
end

end

