data_files = importdata('metadata/data_files.mat'); % list with data folders/files from the base folder
basefolder = 'E:\data\';

% Load data of day and step of each recording
load('metadata/behavior_metadata.mat');

time_in_training = step_data(:,1);
mouse = step_data(:,4);
allSteps = step_data(:,3);
allSteps(allSteps<2) = -1;
allSteps(allSteps==2) = 0;
allSteps(allSteps>2) = 1;

% Calculate the rewarded trials and ratio 
allRatios = [];
allRew = [];
for i = 1:size(data_files, 1)
    fprintf("%d/%d\n", i, size(data_files,1))
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [trimfolder,~] = fileparts(folder);
    % Calculate the control - rewarded licks; % licks in control trials and
    % % licks in reward trials
    [conRewRatio, contLicksPerc,rewLicksPerc] =  getBehaviour(trimfolder);
    allRatios = [allRatios, conRewRatio];
    allRew = [allRew, rewLicksPerc];
end


% Performance over days
data = [allRew',step_data(:,3),step_data(:,2)];

save('data/perf_over_days.mat', 'data')

% Steps over days
data = [step_data(:,3),step_data(:,4),step_data(:,2)];

save('data/steps_over_days.mat', 'data')


%% Performance per environment

perc_correct = nan(size(data_files,1),4);
perc_disp_all= nan(size(data_files,1),2);
for i = 1:size(data_files,1)
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [trimfolder,~] = fileparts(folder);
    [licked, nSal,nNov, displacement] =  getBehaviourpTrav(trimfolder);
    %% correct choice per env
    for s = 1:2
        for n = 1:2
            if n==1
                correct = licked==0;
            else
                correct = licked==1;
            end
            idx = (n-1)*2+s;
            perc_correct(i,idx) = nanmean(correct(nSal == (s-1) & nNov == (n-1)));
        end   
    end
    % Save by displacement
    perc_disp_all(i,1) = nanmean(licked(abs(displacement)<=40 & nNov == 1));
    perc_disp_all(i,2) = nanmean(licked(abs(displacement)>80 & nNov == 1));
end

datatemp = [save_to_R(perc_correct(allSteps<0,:)); save_to_R(perc_correct(allSteps>0,:))];
training = [save_to_R(ones(sum(allSteps<0),4)); save_to_R(ones(sum(allSteps>0),4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4)); save_to_R(repmat(mouse(allSteps>0),1,4))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/perf_per_env.mat', 'data')


datatemp = [save_to_R(perc_disp_all(allSteps<0,:)); save_to_R(perc_disp_all(allSteps>0,:))];
training = [save_to_R(ones(sum(allSteps<0),2)); save_to_R(ones(sum(allSteps>0),2)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,2)); save_to_R(repmat(mouse(allSteps>0),1,2))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/all_perc_p_disp.mat', 'data')


%% Speed per environment
nBins = 100;
vel_ave = nan(size(data_files,1), 4);
vel_ave_p_loc = nan(size(data_files,1), 4,nBins);
for i = 1:size(data_files,1)
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    %% Import and clean loco data
    allLoc = importdata(fullfile(f,'locomotion.mat'));
    inc_loc = allLoc>50/295 & allLoc <= 250/295;
    velRes = getVelocity(allLoc, length(allLoc));
  
    envDown = getEnvVals (f,length(allLoc));
    env = envDown(1,:);
    env = round(env);
    DO = envDown(4,:) ~= 0 & envDown(4,:) ~= 100;
    DO = DO+1; % 1 is no DOel, 2 is DOel
    
    %% Calculate mean velocity
    for n = 1:2 % DOel object or not
        for s = 1:2 % MSC environmnet or not
             idx = (n-1)*2+s;
            incEnv = env==s & DO==n;
            inc = incEnv & inc_loc';
            vel_ave(i,idx) = nanmean(velRes(inc));
            
            % per bin
            loco =   (allLoc(inc));
                loco = loco-min(loco);
                loco = loco./max(loco);
            
                bin_loco = discretize(loco,nBins);
                vel_temp = velRes(inc);
                for loc = 1:nBins
                    vel_ave_p_loc(i,idx,loc) = nanmean(vel_temp(bin_loco==loc));
                end
        end
    end
end

vel_cms = vel_ave*295/7.51; % Calculate back into cm
datatemp = [save_to_R(vel_cms(allSteps<0,1:4)); save_to_R(vel_cms(allSteps>0,1:4))];
training = [save_to_R(ones(sum(allSteps<0),4)); save_to_R(ones(sum(allSteps>0),4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4)); save_to_R(repmat(mouse(allSteps>0),1,4))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/speed_p_env.mat', 'data')

%% Speed relative to parts of the environment

diff_bins = nan(size(data_files,1),4,3);
% familiar location blue ball
diff_bins(:,:,2) = nanmean(vel_ave_p_loc(:,:,16:35),3);
% familiar location control object
diff_bins(:,:,1) = nanmean(vel_ave_p_loc(:,:,66:85),3);
% reward location
diff_bins(:,:,3) = nanmean(vel_ave_p_loc(:,:,86:100),3);
inc_control = [36:65];
diff_bins(:,:,4) = nanmean(vel_ave_p_loc(:,:,inc_control),3);

datatemp = [save_to_R(diff_bins(allSteps<0,1:4,:)); save_to_R(diff_bins(allSteps>0,1:4,:))];
training = [save_to_R(ones(sum(allSteps<0),4,4)); save_to_R(ones(sum(allSteps>0),4,4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4,4)); save_to_R(repmat(mouse(allSteps>0),1,4,4))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/vel_p_bin.mat', 'data')


