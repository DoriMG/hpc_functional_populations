function [pcLoc, trav] = placefields_by_traversal(dataNorm, loc, nLocs, run_frames, trav)
    
    nRois = size(dataNorm,1);
    if size(loc,1) < size(loc,2)
        loc = loc';
    end
    loc = round(loc*nLocs);
    
    % If traversals are not given, than calculate them
    if nargin < 5
          trav = findTraversals(loc, 0.44*max(loc));
    end
    
    ts = max(trav); % Number of traversals
    pcLoc = zeros(nRois, ts,nLocs);

    % Find fluorescence per cell per traversal
    for i = 1:nLocs
        for j = 1:ts
            pcLoc(:,j,i) = nanmean(dataNorm(:, loc==i & trav'==j & run_frames),2);
        end
    end

end