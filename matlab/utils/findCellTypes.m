function [allPC_ps, allOVCs_ps, allPCs, allOVCs] = findCellTypes(f)
%% Variables
track_length = 295;
env_length = 200;
frame_rate = 7.51;
bin_sz = 2;

%% Import and clean 2P data
Fall = importdata(fullfile(f, 'Fall.mat'));
traces = Fall.F-(0.7*Fall.Fneu); % Recommended by Suite2P
traces = traces(logical(Fall.iscell(:,1)),:);
traces = findDeltaF(traces, 1:100);

% Import and preprocess locomotion
allLoc = importdata(fullfile(f,'locomotion.mat'));

% Exclude tunnels
inc_loc = allLoc>50/295 & allLoc <= 250/295;
% Cut to line up data from Arduino and data from 2P acquisition (start time
% synced, but not stop time)
incTime = getIncTime(f, inc_loc);
inc_loc = inc_loc(incTime); 
% normalize location
allLoc = allLoc(incTime);
allLoc = allLoc-min(allLoc);
allLoc = allLoc/max(allLoc);
allLoc = round(allLoc*track_length);

runFrames = getRunFrames(allLoc, track_length, frame_rate); % get running frames from raw locomotion

% get frames for fluorescence
normTraces = traces(:,incTime);


% load environmental data
envVals = getEnvVals (f,length(incTime));
env = round(envVals(1,:));
env = env(incTime); % Sparse (1) vs abundant (2)

% Find object displacement (DO) trials
DO = envVals(4,:) ~= 0 & envVals(4,:) ~= 100; % if this is listed as 0 or 100 the object was not displaced
DO = DO+1; % 1 is no DOel, 2 is DOel
DO = DO(incTime);


% Get location relative to objects
locOfObj = envVals(4,incTime);
locOfObj(locOfObj==0) = 100;
objLoc = allLoc'-(locOfObj);
incObjLoc = objLoc > (-100)& objLoc < (100);


%% find place and OVC cells
allPCs = [];
allOVCs = [];
allPC_ps = [];
allOVCs_ps = [];
allPC_centers = [];
allOVCs_centers = [];
for n = 1:2 % DOel object or not
    for s = 1:2 % sparse or abundant
        inc = env==s & DO==n ;
        inc_pc = inc& inc_loc';
        if sum(inc_pc)>=50
            % Select frames where animal is in the right traversal
            tracesN = normTraces(:,inc_pc);
            loc_temp = allLoc(inc_pc)-50;
            % get running periods
            speed_filtered_frames = runFrames(inc_pc);
            
            
            %% Find place cells
            [PC,p, ~, ~, pfc] = ...
                get_pcs(tracesN, loc_temp, speed_filtered_frames, 7.51, env_length, bin_sz);
        else
            PC = NaN(size(normTraces,1),1);
            p = NaN(size(normTraces,1),1);
            pfc = NaN(size(normTraces,1),1);
        end
        if sum(incObjLoc(inc))>=50
            %% get OVCs
            
            speed_filtered_frames = runFrames(inc);
            tracesN = normTraces(:,inc);
            objLoc_temp = objLoc(inc)'+100;
             
            [OVC,ovc_p, ~, ~, ovc_pfc] = ...
                get_pcs(tracesN(:,incObjLoc(inc)), objLoc_temp(incObjLoc(inc)), speed_filtered_frames(incObjLoc(inc)), 7.51, env_length, bin_sz);
        else
            OVC = NaN(size(normTraces,1),1);
            ovc_p = NaN(size(normTraces,1),1);
            ovc_pfc = NaN(size(normTraces,1),1);
        end
        
        allPC_ps = [allPC_ps, p];
        allPCs = [allPCs, PC];
        allPC_centers = [allPC_centers, pfc];
        
        allOVCs_ps = [allOVCs_ps, ovc_p];
        allOVCs = [allOVCs, OVC];
        allOVCs_centers = [allOVCs_centers, ovc_pfc];
    end
end

save(fullfile(f, 'celltypes.mat'), 'allPC_ps', 'allPCs', 'allPC_centers', 'allOVCs_ps', 'allOVCs', 'allOVCs_centers')
end