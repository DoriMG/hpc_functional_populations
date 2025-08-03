function [pcLoc, pcLocSm] = placefieldsBayes(dataNorm, loc ,nLocs)
%     locdiff = 10;
    nRois = size(dataNorm,1);
    if size(loc,1) < size(loc,2)
        loc = loc';
    end
    
    pcLoc = zeros(nRois, nLocs);
    for i = 1:nLocs
        pcLoc(:,i) = nanmean(dataNorm(:, loc==i),2);
    end
    
    pcLocSm = zeros(size(pcLoc));
    for i = 1:nRois
        pcLocSm(i,:) = smooth(pcLoc(i,:),5);
    end
end