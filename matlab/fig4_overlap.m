
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

% environment data
nBins = 100;
track_length=295;
env_length = 200;
frame_rate=7.51;
bin_sz = env_length/nBins;


%% Fig 4B overlap between conditions
pc_p_folder_overlap= nan(size(data_files,1), 4,4);
pc_p_folder_expected= nan(size(data_files,1), 4,4);
for i = 1:size(data_files,1)
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    %% Import and clean data
    load(fullfile(f,'celltypes.mat'));
    pcs = allPCs;
    
    for idx = 1:4
        overlap = sum(pcs==1&pcs(:,idx)==1)./size(pcs,1);
        pc_p_folder_overlap(i,idx,:)= overlap;
    end
    
    percs = sum(pcs)./size(pcs,1);
    pc_p_folder_expected(i,:,:) = percs.*percs';
end


% Fig 3B - raw overlap
overlap_data  = [nanmean([pc_p_folder_overlap(:,1,3), pc_p_folder_overlap(:,2,4)],2),...
    nanmean([pc_p_folder_overlap(:,1,2), pc_p_folder_overlap(:,3,4)],2),...
    nanmean([pc_p_folder_overlap(:,1,4), pc_p_folder_overlap(:,2,3)],2)];
datatemp = [save_to_R(overlap_data(allSteps<0,:,:)); save_to_R(overlap_data(allSteps>0,:,:))];
training = [save_to_R(ones(sum(allSteps<0),3,1)); save_to_R(ones(sum(allSteps>0),3,1)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,3,1)); save_to_R(repmat(mouse(allSteps>0),1,3,1))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/pc_overlap_raw.mat', 'data')

% Fig 3B - Overlap compared to chance
new_overlap  = (pc_p_folder_overlap-pc_p_folder_expected)./pc_p_folder_expected;
overlap_data  = [nanmean([new_overlap(:,1,3), new_overlap(:,2,4)],2),...
    nanmean([new_overlap(:,1,2), new_overlap(:,3,4)],2),...
    nanmean([new_overlap(:,1,4), new_overlap(:,2,3)],2)];
datatemp = [save_to_R(overlap_data(allSteps<0,:,:)); save_to_R(overlap_data(allSteps>0,:,:))];
training = [save_to_R(ones(sum(allSteps<0),3,1)); save_to_R(ones(sum(allSteps>0),3,1)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,3,1)); save_to_R(repmat(mouse(allSteps>0),1,3,1))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/pc_overlap_p_cond_diff.mat', 'data')



%% Fig 4F-G - overlap in place fields
nBins = 100;

allAveOverlap = nan(size(data_files,1), 4,4);
allAveOverlap_pc = nan(size(data_files,1), 4,4);
for i = 1:size(data_files,1)
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    load(fullfile(f,'celltypes.mat'));
    pcs = allPCs;
    load(fullfile(f, 'data.mat'));
    
    spks = traces;
    %% First extract all place fields
    all_pcLoc = nan(size(spks,1),nBins, 4);
    for n = 1:2 % DOel object or not
        for s = 1:2 % MSC environmnet or not
            idx = (n-1)*2+s;
            incEnv = env==s & DO==n;
            inc = incEnv'&inc_loc&runFrames;
            loco = (allLoc(inc) - min(allLoc(inc)))/(max(allLoc(inc))- min(allLoc(inc)));
            
            [pcLoc, tune_all]= placefieldsBayes(spks(:,inc), round(loco*nBins), nBins);
            all_pcLoc(:,:,idx) = pcLoc;
        end
    end
    % Now repeat and correlate between all combinations of conditions
    for idx1 = 1:4
        for idx2 = 1:4
            corr_temp = [];
            corr_pc_temp = [];
            for c = 1:size(pcLoc,1)
                inc = ~isnan(all_pcLoc(c,:,idx1))&~isnan( all_pcLoc(c,:,idx2));
                if sum(inc)>5 % Have at least 5 known bins
                    cc = corrcoef(all_pcLoc(c,inc,idx1), all_pcLoc(c,inc,idx2));
                    corr_temp = [corr_temp, cc(1,2)];
                    % This is only PCs
                    if pcs(c,idx1)==1 && pcs(c,idx2)==1
                        corr_pc_temp = [corr_pc_temp, cc(1,2)];
                    end
                end
                
            end
            allAveOverlap(i,idx1,idx2) = nanmean(corr_temp);
            allAveOverlap_pc(i,idx1,idx2) = nanmean(corr_pc_temp);
        end
    end
    
end
% Fig 4F - all cells
overlap_data  = cat(2,cat(1,allAveOverlap(:,1,3), allAveOverlap(:,2,4)),...
    cat(1,allAveOverlap(:,1,2), allAveOverlap(:,3,4)),...
    cat(1,allAveOverlap(:,1,4), allAveOverlap(:,2,3)));
datatemp = save_to_R(overlap_data);
training = save_to_R(repmat(allSteps, 2, 3));
mice =save_to_R(repmat(mouse, 2, 3));
data = [datatemp, training(:,1), mice(:,1)];
save('data/all_overlap_corr.mat', 'data')

% Fig 4G - PCs
overlap_data  = cat(2,cat(1,allAveOverlap_pc(:,1,3), allAveOverlap_pc(:,2,4)),...
    cat(1,allAveOverlap_pc(:,1,2), allAveOverlap_pc(:,3,4)),...
    cat(1,allAveOverlap_pc(:,1,4), allAveOverlap_pc(:,2,3)));
datatemp = save_to_R(overlap_data);
training = save_to_R(repmat(allSteps, 2, 3));
mice =save_to_R(repmat(mouse, 2, 3));
data = [datatemp, training(:,1), mice(:,1)];
save('data/pc_overlap_corr.mat', 'data')

%% Run decoder
rng(10) % Set for reproducibility
nbins = 100;

all_error_p_env_split = nan(size(data_files,1),4, 4);
for i = 1:size(data_files,1)
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    load(fullfile(f, 'data.mat'));
    load(fullfile(f,'celltypes.mat'));
    allPCs(isnan(allPCs)) = 0;
    
    load(fullfile(f,'Fall.mat'));
    spks = spks(logical(iscell(:,1)),:)>0;

    tras = findTraversals(allLoc, 0.5*max(allLoc));
    for n = 1:4 % Run through all conditions
        for s = 1:4 % Run through all conditions for comparison
            % get conditions from index
            [s_env, do_n] = ind2sub([2 2], n);
            incEnv = env==s_env & DO==do_n;

            inc = incEnv'&inc_loc&runFrames;
                
            if n == s % Decode from self by taking test and training set, split by traversals
                    tras_list = unique(tras(inc));
                    rand_tras = tras_list(randperm(length(tras_list)));
                    inc_tras = ismember(tras,rand_tras(1:ceil(length(rand_tras)/2)));
                    inc_tras_test = ismember(tras,rand_tras(ceil(length(rand_tras)/2)+1:end));
            else % If decode from other environment, can include all traversals
                inc_tras = ones(size(inc_tras));
                inc_tras_test = ones(size(inc_tras));
            end
            
            inc_train =  incEnv'&inc_loc&runFrames & inc_tras';
            inc_test =  incEnv'&inc_loc&runFrames & inc_tras_test';
            loco =   (allLoc(inc_train));
            loco_test = (allLoc(inc_test));
            
            [pcLoc, tune_all]= placefieldsBayes(spks(:,inc_train), round(loco*nbins), nbins);
            pAgs = nan(nbins, 1, size(tune_all,1));
            for b = 1:size(tune_all,1)
                pAgs(:,:,b) = tune_all(b,:);
            end
            
            % Now get the data from comparison environment and decode
            [s_env, do_n] = ind2sub([2 2], s);
            incEnv = env==s_env & DO==do_n;
            if n == s
               inc = inc_test; 
            else
                inc =  incEnv'&inc_loc&runFrames;
            end
            loco =   (allLoc(inc));
            [ post ] = decode_calcPost_twoPhoton( spks(:,inc),  pAgs);
            probs = squeeze(post)';
            true_loco = discretize(loco,nbins);
            [~,pred] = max(probs,[],2);
            error = min(abs(pred-true_loco), 147.5-abs(pred-true_loco));

            all_error_p_env_split(i,n,s) = nanmean(error);
           
        end
    end
end

overlap_data  = cat(2,cat(1,all_error_p_env_split(:,1,3), all_error_p_env_split(:,2,4),all_error_p_env_split(:,3,1), all_error_p_env_split(:,4,2)),...
    cat(1,all_error_p_env_split(:,1,2), all_error_p_env_split(:,3,4),all_error_p_env_split(:,2,1), all_error_p_env_split(:,4,3)),...
    cat(1,all_error_p_env_split(:,1,4), all_error_p_env_split(:,2,3),all_error_p_env_split(:,4,1), all_error_p_env_split(:,3,2)));

datatemp = save_to_R(overlap_data);
training = save_to_R(repmat(allSteps, 4, 3));
mice =save_to_R(repmat(mouse, 4, 3));
data = [datatemp, training(:,1), mice(:,1)];
save('data/all_cross_decoder_error.mat', 'data')

% Subtract the control (i.e. environment decoded from itself
overlap_data  = cat(2,cat(1,all_error_p_env_split(:,1,3) - all_error_p_env_split(:,3,3), ...
    all_error_p_env_split(:,2,4)- all_error_p_env_split(:,4,4), ...
    all_error_p_env_split(:,3,1)- all_error_p_env_split(:,1,1), ...
    all_error_p_env_split(:,4,2)- all_error_p_env_split(:,2,2)),...
    cat(1,all_error_p_env_split(:,1,2)- all_error_p_env_split(:,2,2), ...
    all_error_p_env_split(:,3,4)- all_error_p_env_split(:,4,4), ...
    all_error_p_env_split(:,2,1)- all_error_p_env_split(:,1,1), ...
    all_error_p_env_split(:,4,3)- all_error_p_env_split(:,3,3)),...
    cat(1,all_error_p_env_split(:,1,4)- all_error_p_env_split(:,4,4), ...
    all_error_p_env_split(:,2,3)- all_error_p_env_split(:,3,3), ...
    all_error_p_env_split(:,4,1)- all_error_p_env_split(:,1,1), ...
    all_error_p_env_split(:,3,2)- all_error_p_env_split(:,2,2)));

datatemp = save_to_R(overlap_data);
training = save_to_R(repmat(allSteps, 4, 3));
mice =save_to_R(repmat(mouse, 4, 3));
data = [datatemp, training(:,1), mice(:,1)];
save('data/all_cross_decoder_diff.mat', 'data')

