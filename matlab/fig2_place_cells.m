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

% Fig 2B number of place cells
pc_p_folder = nan(size(data_files,1), 4);
for i = 1:size(data_files,1)
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    %% Import and clean data
    load(fullfile(f,'celltypes.mat'));
    pcs = allPCs;
    
    pc_p_folder(i,:) = sum(pcs)/size(pcs,1);
end

datatemp = [save_to_R(pc_p_folder(allSteps<0,:)); save_to_R(pc_p_folder(allSteps>0,:))];
training = [save_to_R(ones(sum(allSteps<0),4)); save_to_R(ones(sum(allSteps>0),4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4)); save_to_R(repmat(mouse(allSteps>0),1,4))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/mean_pc_perc.mat', 'data')



%% Fig 2C Location of place cells

nBins = 20;
pc_center_per_loc = nan(size(data_files,1), 4,nBins);
pc_center_per_loc_nopc = nan(size(data_files,1), 4,nBins);
for i = 1:size(data_files,1)
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    load(fullfile(f, 'data.mat'));
    load(fullfile(f,'celltypes.mat'));
    allPCs(isnan(allPCs)) = 0;
    
    %% Do Bayesian decoder
    for n = 1:2 % DOel object or not
        for s = 1:2 % MSC environmnet or not
            idx = (n-1)*2+s;
                incEnv = env==s & DO==n;
                inc = incEnv'&inc_loc&runFrames;
                tempLoc = (allLoc(inc) - min(allLoc(inc)))/(max(allLoc(inc))- min(allLoc(inc)));
                tempLoc = discretize(tempLoc, nBins);
                tempFluor = traces(:, inc);
                
                [pcLoc, pcLocSm]= placefieldsBayes(tempFluor, tempLoc, nBins);
                [~,pcCentre] = nanmax(pcLoc,[],2);
                for loc = 1:nBins
                    pc_center_per_loc(i,idx, loc) = sum(pcCentre(logical(allPCs(:,idx)))==loc)/sum(logical(allPCs(:,idx)));
                    pc_center_per_loc_nopc(i,idx, loc) = sum(pcCentre(logical(~allPCs(:,idx)))==loc)/size(tempFluor,1);
                end
        end
    end
end

diff_bins = nan(size(data_files,1),4,3);
% familiar location blue ball
diff_bins(:,:,2) = nanmean(pc_center_per_loc(:,:,4:7),3);
% familiar location control object
diff_bins(:,:,1) = nanmean(pc_center_per_loc(:,:,14:17),3);
% reward location
diff_bins(:,:,3) = nanmean(pc_center_per_loc(:,:,18:20),3);
inc_control = [8:13];
diff_bins(:,:,4) = nanmean(pc_center_per_loc(:,:,inc_control),3);

datatemp = [save_to_R(diff_bins(allSteps<0,1:4,:)); save_to_R(diff_bins(allSteps>0,1:4,:))];
training = [save_to_R(ones(sum(allSteps<0),4,4)); save_to_R(ones(sum(allSteps>0),4,4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4,4)); save_to_R(repmat(mouse(allSteps>0),1,4,4))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/pf_centres_binned_pc.mat', 'data')


%% Fig 2 - S Fig 2 PC characteristics

allAvePfs = nan(size(data_files,1), 4);
allAveStab = nan(size(data_files,1), 4);
allAveMI = nan(size(data_files,1), 4);
allIO = nan(size(data_files,1), 4);
for i = 1:size(data_files,1)
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    %% Import and clean data
    load(fullfile(f,'celltypes.mat'));
    pcs = allPCs;
    
    load(fullfile(f, 'data.mat'));

    
    % Loop over conditions
    for n = 1:2 % DOel object or not
        for s = 1:2 % MSC environmnet or not
            idx = (n-1)*2+s;
            incEnv = env==s & DO==n;
            inc = incEnv'&inc_loc&runFrames;
            spks = traces(:,inc);
            
            tempLoc = allLoc(inc);
            tempLoc = tempLoc-min(tempLoc);
            tempLoc = tempLoc/max(tempLoc);
            
            [pcLoc, pcLocSm]= placefieldsBayes(spks, round(tempLoc*nBins), nBins);
         
            [obj_loc, obj_trav]= placefields_by_traversal(spks, tempLoc, 50, runFrames(inc));

            % Calculate MI
            if ~isempty(tempLoc) & ~isnan(tempLoc)
                all_mi = findMutualInfo(spks, tempLoc*track_length, track_length, bin_sz);
            else
                all_mi = nan(size(spks,1),1);
            end
            
            
            if sum(isnan(pcs(:,idx)))==0
                pcs_temp = logical(pcs(:,idx));
                
                % place field size
                [sizePf] = findPlaceFieldSizes(pcLoc(pcs_temp,:), env_length/nBins);
                avePf = nanmean(sizePf);
                
                % stability
                stabs = findStability(obj_loc);
                aveStab = nanmean(stabs(pcs_temp));
                
                % Mutual information
                aveMI = nanmean(all_mi(pcs_temp));
                
                % I/O rate
                [sizePf] = findOutInFieldRate(pcLoc(pcs_temp,:));
                aveIO = nanmean(sizePf);
            else
                avePf = nan;
                aveStab = nan;
                aveMI = nan;
                aveIO = nan;
            end
            allAvePfs(i,idx) = avePf;
            allAveStab(i,idx) = aveStab;
            allAveMI(i,idx) = aveMI;
            allIO(i,idx) = aveIO;

        end
    end
end

datatemp = [save_to_R(allIO(allSteps<0,:)); save_to_R(allIO(allSteps>0,:))];
training = [save_to_R(ones(sum(allSteps<0),4)); save_to_R(ones(sum(allSteps>0),4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4)); save_to_R(repmat(mouse(allSteps>0),1,4))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/allIO.mat', 'data')

datatemp = [save_to_R(allAveMI(allSteps<0,:)); save_to_R(allAveMI(allSteps>0,:))];
training = [save_to_R(ones(sum(allSteps<0),4)); save_to_R(ones(sum(allSteps>0),4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4)); save_to_R(repmat(mouse(allSteps>0),1,4))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/allAveMI.mat', 'data')

datatemp = [save_to_R(allAveStab(allSteps<0,:)); save_to_R(allAveStab(allSteps>0,:))];
training = [save_to_R(ones(sum(allSteps<0),4)); save_to_R(ones(sum(allSteps>0),4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4)); save_to_R(repmat(mouse(allSteps>0),1,4))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/allAveStab.mat', 'data')

datatemp = [save_to_R(allAvePfs(allSteps<0,:)); save_to_R(allAvePfs(allSteps>0,:))];
training = [save_to_R(ones(sum(allSteps<0),4)); save_to_R(ones(sum(allSteps>0),4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4)); save_to_R(repmat(mouse(allSteps>0),1,4))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/allAvePfs.mat', 'data')


%% Fig 2F - S Fig2D -  Average mean error
% environment data
nbins=100;
nReps = 100;
% Run decoder across all cells
all_error_p_env_split = nan(size(data_files,1),4);
all_error_p_env_split_rand = nan(size(data_files,1),4, nReps);
for i = 1:size(data_files,1)
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    load(fullfile(f, 'data.mat'));
    load(fullfile(f,'celltypes.mat'));
    allPCs(isnan(allPCs)) = 0;
    
    load(fullfile(f,'Fall.mat'));
    spks = spks>0;
    
    mean_error_p_env = nan(4,1);
    all_post = {};
    for n = 1:2 % DOel object or not
        for s = 1:2 % MSC environmnet or not
            idx = (n-1)*2+s;
            incEnv = env==s & DO==n;
            inc = incEnv'&inc_loc&runFrames;
            loco = (allLoc(inc) - min(allLoc(inc)))/(max(allLoc(inc))- min(allLoc(inc)));
            
            [pcLoc, tune_all]= placefieldsBayes(spks(:,inc), round(loco*nbins), nbins);
            pAgs = nan(nbins, 1, size(spks,1));
            for b = 1:size(spks,1)
                pAgs(:,:,b) = tune_all(b,:);
            end
            
            [ post ] = decode_calcPost_twoPhoton( spks(:,inc),  pAgs);
            probs = squeeze(post)';
            true_loco = discretize(loco,nbins);
            [~,pred] = max(probs,[],2);
            error = min(abs(pred-true_loco), 149.5-abs(pred-true_loco)); % to account for circular environments
            
            all_error_p_env_split(i,idx) = nanmean(error);
            
            % perform random control by randomizing cell order
             for rep = 1:nReps %
                cell_order = randperm(size(spks,1));
                mean_error_p_env_rand = nan(4,1);
                
                % Decode with randomized cell order
                [ post ] = decode_calcPost_twoPhoton( spks(cell_order,inc),  pAgs);
                probs = squeeze(post)';
                true_loco = discretize(loco,nbins);
                [~,pred] = max(probs,[],2);
                error = min(abs(pred-true_loco), 147.5-abs(pred-true_loco));
                bin_loco = discretize(loco,nLocs);
               
                all_error_p_env_split_rand(i,idx,rep) = nanmean(error);
             end
        end
    end
end

all_error_p_env = all_error_p_env_split .* [2,2,2,2]; % convert to cm
datatemp = [save_to_R(all_error_p_env(allSteps<0,:)); save_to_R(all_error_p_env(allSteps>0,:))];
training = [save_to_R(ones(sum(allSteps<0),4)); save_to_R(ones(sum(allSteps>0),4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4)); save_to_R(repmat(mouse(allSteps>0),1,4))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/all_error.mat', 'data')


all_error_p_env_rand = squeeze(nanmean(all_error_p_env_split_rand,3)) .* [2,2,2,2]; % convert to cm
datatemp = [save_to_R(all_error_p_env_rand(allSteps<0,:)); save_to_R(all_error_p_env_rand(allSteps>0,:))];
training = [save_to_R(ones(sum(allSteps<0),4)); save_to_R(ones(sum(allSteps>0),4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4)); save_to_R(repmat(mouse(allSteps>0),1,4))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/all_error_rand.mat', 'data')

%% Fig 2G - Error pc - rand
pc_error_p_env_split = nan(size(data_files,1),4);
pc_error_p_env_rand = nan(size(data_files,1),4);


for i = 1:size(data_files,1)
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    load(fullfile(f, 'data.mat'));
    load(fullfile(f,'celltypes.mat'));
    allPCs(isnan(allPCs)) = 0;
    
    load(fullfile(f,'Fall.mat'));
    spks = spks>0;
    
    mean_error_p_env = nan(4,1);
    min_PCs = min(sum(allPCs));
    %% Do Bayesian decoder
    if min_PCs>0
        for n = 1:2 % DOel object or not
            for s = 1:2 % MSC environmnet or not
                idx = (n-1)*2+s;
                incEnv = env==s & DO==n;
                inc = incEnv'&inc_loc&runFrames;
                loco = (allLoc(inc) - min(allLoc(inc)))/(max(allLoc(inc))- min(allLoc(inc)));
                
                errors = nan(10,1);
                errors_random = nan(10,1);
                for nrep = 1:10
                    % Actual error
                    pc_idx = find(allPCs(:,idx));
                    pc_inc = randsample(pc_idx,min_PCs);
                    [pcLoc, tune_all]= placefieldsBayes(spks(pc_inc,inc), round(loco*nbins), nbins);
                    pAgs = nan(nbins, 1, size(tune_all,1));
                    for b = 1:size(tune_all,1)
                        pAgs(:,:,b) = tune_all(b,:);
                    end
                    
                    [ post ] = decode_calcPost_twoPhoton( spks(pc_inc,inc),  pAgs);
                    all_post{idx} = post;
                    probs = squeeze(post)';
                    true_loco = discretize(loco,nbins);
                    [~,pred] = max(probs,[],2);
                    error = min(abs(pred-true_loco), 147.5-abs(pred-true_loco));
                    errors(nrep) = nanmean(error);
                    
                    % shuffled error
                    pc_idx = find(~allPCs(:,idx));
                    pc_inc = randsample(pc_idx,min_PCs);
                    [pcLoc, tune_all]= placefieldsBayes(spks(pc_inc,inc), round(loco*nbins), nbins);
                    pAgs = nan(nbins, 1, size(tune_all,1));
                    for b = 1:size(tune_all,1)
                        pAgs(:,:,b) = tune_all(b,:);
                    end
                    
                    [ post ] = decode_calcPost_twoPhoton(spks( pc_inc,inc),  pAgs);
                    all_post{idx} = post;
                    probs = squeeze(post)';
                    true_loco = discretize(loco,nbins);
                    [~,pred] = max(probs,[],2);
                    error = min(abs(pred-true_loco), 147.5-abs(pred-true_loco));
                    errors_random(nrep) = nanmean(error);
                end
                
                pc_error_p_env_split(i,idx) = nanmean(errors);
                pc_error_p_env_rand(i,idx) = nanmean(errors_random);
                
            end
        end
    end
    
end

all_error_p_env = pc_error_p_env_split .* [2,2,2,2]; % convert to cm
all_error_p_env_rand = pc_error_p_env_rand .* [2,2,2,2]; % convert to cm
datatemp = [save_to_R(all_error_p_env(allSteps<0,:)); save_to_R(all_error_p_env(allSteps>0,:))];
rand_data = [save_to_R(all_error_p_env_rand(allSteps<0,:)); save_to_R(all_error_p_env_rand(allSteps>0,:))];
training = [save_to_R(ones(sum(allSteps<0),4)); save_to_R(ones(sum(allSteps>0),4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4)); save_to_R(repmat(mouse(allSteps>0),1,4))];
data = [datatemp, training(:,1), mice(:,1)];
data(:,1) = data(:,1)-rand_data(:,1);
save('data/pc_error.mat', 'data')
