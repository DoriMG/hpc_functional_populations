data_files = importdata('metadata/data_files.mat'); % list with data folders/files from the base folder
basefolder = 'D:\Frankfurt\Papers\Sussex\hippocampal_pops\v3\data';

% Load data of day and step of each recording
load('metadata/behavior_metadata.mat');

time_in_training = step_data(:,1);
mouse = step_data(:,4);
allSteps = step_data(:,3);
allSteps(allSteps<2) = -1;
allSteps(allSteps==2) = 0;
allSteps(allSteps>2) = 1;

track_length = 295;
env_length = 200;
frame_rate = 7.51;


for i = 1:size(data_files,1)
    fprintf("%d/%d\n", i, size(data_files,1))
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    %% Import and clean fluorescence data
    Fall = importdata(fullfile(f, 'Fall.mat'));
    traces = Fall.F-(0.7*Fall.Fneu);
    traces = traces(logical(Fall.iscell(:,1)),:);
    traces = findDeltaF(traces, 1:100);
    
    % Cut to line up data from Arduino and data from 2P acquisition (start time
    % synced, but not stop time)
    incTime = getIncTime(f, inc_loc);
    
     %% Import and clean locomotion fluorescence data
    allLoc = importdata(fullfile(f,'locomotion.mat'));
    inc_loc = allLoc>50/295 & allLoc <= 250/295;
    inc_loc = inc_loc(incTime);
    
    % normalize location
    allLoc = allLoc(incTime);
    allLoc = round(allLoc*track_length); % convert to cm

    runFrames = getRunFrames(allLoc, track_length, frame_rate); % get running frames from raw locomotion

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
    
    % Save all in one big mat file
    
    save(fullfile(f, 'data.mat'), 'allLoc', 'traces', 'env', 'DO', 'objLoc', 'incObjLoc', 'runFrames', 'inc_loc')
end