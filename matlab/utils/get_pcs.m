function [PC,p, rms, maxInt, pfc] = get_pcs(traces, loc, speed_filtered_frames, frame_rate, track_length, bin_sz)
%inputs 
% traces     -    nRois x frames array containing the delta F/F fluorescence
%                 traces
% loc        -    1 x frames location trace (cm) 
% speed_filtered_frames 
%            -    1 x frames binary array where 1s are frames with speed
%                 high enough to retain the 
% frame_rate -    double, frame rate of recording (fps)
% track_length
%            -    int, length of the track (cm)
% bin_sz     -    int, size of the bin (cm)

nReps = 500; % Number of shuffles

nROIs = size(traces,1);

% Calculate the ratemaps for each cell
[~,binnedPos,binsUnique]=histnd2(loc(speed_filtered_frames), track_length, bin_sz); 
df_f = traces';
rms = nan(nROIs, ceil(track_length/bin_sz));
for ROI=1:nROIs
    [ratemap,~,~,~,~,~]= find_fields(ROI,binnedPos, binsUnique,df_f(speed_filtered_frames,ROI),1,bin_sz);  
    rms(ROI,:) = ratemap;
end

% Calculate the peak value for each cell
[maxInt, pfc] = max(rms,[],2);


maxValues = nan(size(traces,1), nReps); % Will hold peak values for shuffles

% Circularly shift traces in time 500 times by a random period > 5 s
for reps = 1:nReps
    % Get random offset (> 5 seconds) (and make sure it's not more than the
    % total length - 5 seconds)
    shufMin = min(ceil(frame_rate*5), round(size(traces,2)/4));
    shufMax = size(traces,2)-2*shufMin;
    randOffset = randi(shufMax)+shufMin;
    
    % Shuffle traces and calculate rms
    shuffTrace = circshift(traces,randOffset,2);
    
    % Recalculate ratemaps and get maximum value
    df_f = shuffTrace';
    for ROI=1:nROIs
        [ratemap,rmErr]= find_rms(binnedPos,binsUnique, df_f(speed_filtered_frames,ROI),bin_sz);
        rmsshuff(ROI,:) = ratemap';
    end
    % Find the peak value for shuffle
    maxValues(:,reps) = max(rmsshuff,[],2);
end

% Find the 99 percentile threshold for the peak value for each cell
thresh = nan(nROIs,1);
p = nan(nROIs,1);
for roi = 1:nROIs
    % Calculate 99 percentile threshold
    thresh(roi) = prctile(maxValues(roi,:),99);

    try
        % Calculates the percentile value of the peak compared to the shuffles
        p(roi) = invprctile(maxValues(roi,:), maxInt(roi))/100;
    catch
        p(roi) = NaN;
    end
end
totInt = nansum(traces, 2); % Find the total activity of each cell

% Make sure cells with no activity are not selected as PCs
p(totInt == 0) = 0;
thresh(totInt == 0) = Inf;

% Cell is PC if the peak value is higher than the threshold
PC =  maxInt>=thresh;

