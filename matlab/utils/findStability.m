function stability = findStability(pcLoc)

% Define even and odd traversals
nTravs = size(pcLoc,2);
evenHalf = 2:2:nTravs;
oddHalf = 1:2:nTravs;

% Find mean  fluorescence in even vs odd traversals
avePcLocEven = squeeze(nanmean(pcLoc(:,evenHalf,:),2));
avePcLocOdd = squeeze(nanmean(pcLoc(:,oddHalf,:),2));

% Calculate the stability per cell
stability = nan(size(pcLoc,1),1);
for cell = 1:size(pcLoc,1)
    exc = isnan(avePcLocEven(cell,:)) | isnan(avePcLocOdd(cell,:)); % Make sure there's enough locations with soem fluorescence
    if sum(exc) < size(pcLoc,3)/2
        cc = corrcoef(avePcLocEven(cell,~exc), avePcLocOdd(cell,~exc));
        stability(cell) = cc(1,2);
    else
        stability(cell) = nan;
    end
end