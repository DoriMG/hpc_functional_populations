% Note: For ease everything is calculated over 4 environments, but since
% they're only formally defined in shifted environments, so the familiar
% will not be used for plotting/stats

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

% Fig 3B number of OVCs
pc_p_folder = nan(size(data_files,1), 4);
for i = 1:size(data_files,1)
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    %% Import and clean data
    load(fullfile(f,'celltypes.mat'));
    allPCs(isnan(allPCs)) = 0;
    allOVCs(isnan(allOVCs)) = 0;
    pcs = logical(allPCs);
    ovc = logical(allOVCs);
    
    pc_p_folder(i,:) = sum(ovc & pcs)/size(pcs,1);
    
end

datatemp = [save_to_R(pc_p_folder(allSteps<0,:)); save_to_R(pc_p_folder(allSteps>0,:))];
training = [save_to_R(ones(sum(allSteps<0),4)); save_to_R(ones(sum(allSteps>0),4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4)); save_to_R(repmat(mouse(allSteps>0),1,4))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/mean_ovc_perc.mat', 'data')



%% Fig 3 - S Fig 4 OVC characteristics

allAvePfs = nan(size(data_files,1), 4);
allAveStab = nan(size(data_files,1), 4);
allAveMI = nan(size(data_files,1), 4);
allIO = nan(size(data_files,1), 4);
for i = 1:size(data_files,1)
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    %% Import and clean data
    load(fullfile(f,'celltypes.mat'));
    allPCs(isnan(allPCs)) = 0;
    allOVCs(isnan(allOVCs)) = 0;
    pcs = logical(allPCs);
    ovc = logical(allOVCs);
    pcs = ovc & ~pcs; % OVCs is any cells that are ovcs, but not pcs
    
    load(fullfile(f, 'data.mat'));
    
    allLoc = objLoc'; % Use the location relative to the object
    % Loop over conditions
    for n = 1:2 % DOel object or not
        for s = 1:2 % MSC environmnet or not
            idx = (n-1)*2+s;
            incEnv = env==s & DO==n;
            inc = incEnv'&runFrames&incObjLoc';
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
save('data/allIO_ovc.mat', 'data')

datatemp = [save_to_R(allAveMI(allSteps<0,:)); save_to_R(allAveMI(allSteps>0,:))];
training = [save_to_R(ones(sum(allSteps<0),4)); save_to_R(ones(sum(allSteps>0),4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4)); save_to_R(repmat(mouse(allSteps>0),1,4))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/allAveMI_ovc.mat', 'data')

datatemp = [save_to_R(allAveStab(allSteps<0,:)); save_to_R(allAveStab(allSteps>0,:))];
training = [save_to_R(ones(sum(allSteps<0),4)); save_to_R(ones(sum(allSteps>0),4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4)); save_to_R(repmat(mouse(allSteps>0),1,4))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/allAveStab_ovc.mat', 'data')

datatemp = [save_to_R(allAvePfs(allSteps<0,:)); save_to_R(allAvePfs(allSteps>0,:))];
training = [save_to_R(ones(sum(allSteps<0),4)); save_to_R(ones(sum(allSteps>0),4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4)); save_to_R(repmat(mouse(allSteps>0),1,4))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/allAvePfs_ovc.mat', 'data')



%% Figure 3 F - OVC locations
nBins = 20;

pc_p_folder = nan(size(data_files,1), 4);
pf_centres = nan(size(data_files,1),4,nBins);
mean_fluor = nan(size(data_files,1),4,nBins);
for i = 1:size(data_files,1)
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    
    load(fullfile(f,'celltypes.mat'));
    allPCs(isnan(allPCs)) = 0;
    allOVCs(isnan(allOVCs)) = 0;
    pcs = logical(allPCs);
    ovc = logical(allOVCs);
    pcs = ovc & ~pcs; % OVCs is any cells that are ovcs, but not pcs
    
    load(fullfile(f, 'data.mat'));

    for n = 1:2 % DOel object or not
        for s = 1:2 % MSC environmnet or not
            idx = (n-1)*2+s;
            incEnv = env==s & DO==n;
            inc = incEnv'&inc_loc&runFrames&incObjLoc';
            spks = traces(:,inc);
            ovc_cells = pcs(:,idx);
            
            tempLoc = allLoc(inc);
            tempLoc = tempLoc-min(tempLoc);
            tempLoc = tempLoc/max(tempLoc);
            
            [pcLoc, pcLocSm]= placefieldsBayes(spks, round(tempLoc*nBins), nBins);
            
            mean_fluor(i,idx,:) = nanmean(pcLoc(ovc_cells,:),1);
            [~,centre_temp] = nanmax(pcLoc(ovc_cells,:),[],2);
            for nb = 1:nBins
                pf_centres(i,idx,nb) = nansum(centre_temp==nb)./sum(ovc_cells);
            end
        end
    end
end

diff_bins = nan(size(data_files,1),4,3);
% familiar location blue ball
diff_bins(:,:,1) = nansum(pf_centres(:,:,1:7),3)/7/10;
% familiar location control object
diff_bins(:,:,2) = nansum(pf_centres(:,:,8:13),3)/6/10;
% reward location
diff_bins(:,:,3) = nansum(pf_centres(:,:,14:20),3)/7/10;
datatemp = [save_to_R(diff_bins(allSteps<0,1:4,:)); save_to_R(diff_bins(allSteps>0,1:4,:))];
training = [save_to_R(ones(sum(allSteps<0),4,3)); save_to_R(ones(sum(allSteps>0),4,3)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4,3)); save_to_R(repmat(mouse(allSteps>0),1,4,3))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/pf_centres_binned_ovc.mat', 'data')


%% Fig 3E Correlation
nBins = 20;
nReps = 100;
corr_p_env = nan(size(data_files,1), 4);
corr_p_env_rand= nan(size(data_files,1), 4);
for i = 1:size(data_files,1)
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    %% Import and clean data
    load(fullfile(f,'celltypes.mat'));
    allPCs(isnan(allPCs)) = 0;
    allOVCs(isnan(allOVCs)) = 0;
    pcs = logical(allPCs);
    ovc = logical(allOVCs);
    pcs = ovc & ~pcs; % OVCs is any cells that are ovcs, but not pcs
    
    load(fullfile(f, 'data.mat'));
    tras = findTraversals(allLoc, 0.5*max(allLoc));
    
    
    inc_loc = allLoc>100 & allLoc<300; % 50 around control obj
    
    for n = 2 % DOel object or not
        for s = 1:2 % MSC environmnet or not
            idx = (n-1)*2+s;
            incEnv = env==s & DO==n;
            
            % Map around cue object
            inc = incEnv'&runFrames&incObjLoc';
            spks = traces(:,inc);
            ovc_cells = pcs(:,idx);
            
            tempLoc = objLoc(inc);
            tempLoc = tempLoc-min(tempLoc);
            tempLoc = tempLoc/max(tempLoc);
            
            [pcLocObj, pcLocSm]= placefieldsBayes(spks, round(tempLoc*nBins), nBins);
            
            % Map around control ojbect
            inc = incEnv'&runFrames & inc_loc;
            spks = traces(:,inc);
            tempLoc = allLoc(inc)';
            tempLoc = tempLoc-min(tempLoc);
            tempLoc = tempLoc/200;
            
            [pcLoc, pcLocSm]= placefieldsBayes(spks, round(tempLoc*nBins), nBins);
            
            ovc_cells = find(ovc_cells);
            corr_temp = [];
            corr_temp_rand = [];
            corr_all_cells = nan(length(ovc_cells),1);
            if length(ovc_cells)>1 % Need at least 2 cells to d othe shuffle
                for c = 1:length(ovc_cells)
                    incTemp = ~isnan(pcLocObj(ovc_cells(c),:)) & ~isnan(pcLoc(ovc_cells(c),:));
                    if sum(incTemp)>15 % have at least 15 valid frames ( 2 secs)
                        cc = corrcoef(pcLocObj(ovc_cells(c),incTemp), pcLoc(ovc_cells(c),incTemp));
                        corr_temp = [corr_temp, cc(1,2)];
                    else
                        corr_temp = NaN;
                    end
                    
                    
                    for rep = 1:nReps
                        rno = randi(size(pcLoc,1)); % select random cells
                        while rno == c % Do not want to include self comparisons in the shuffle
                            rno = randi(size(pcLoc,1));
                        end
                        incTemp = ~isnan(pcLocObj(ovc_cells(c),:)) & ~isnan(pcLoc(rno,:));
                        
                        
                        if sum(incTemp)>15
                            cc = corrcoef(pcLocObj(ovc_cells(c),incTemp), pcLoc(rno,incTemp));
                            corr_temp_rand = [corr_temp_rand, cc(1,2)];
                        else
                            corr_temp_rand = NaN;
                        end
                    end
                    corr_all_cells(c) = nanmean(corr_temp_rand);
                end
                corr_p_env(i,idx) = nanmean(corr_temp);
                corr_p_env_rand(i,idx) = nanmean(corr_all_cells);
            end
        end
    end
end


temp = cat(3,corr_p_env(:,3:4), corr_p_env_rand(:,3:4)); % Only include shifted environments
datatemp = [save_to_R(temp(allSteps<0,:,:)); save_to_R(temp(allSteps>0,:,:))];
training = [save_to_R(ones(sum(allSteps<0),2,2)); save_to_R(ones(sum(allSteps>0),2,2)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,2,2)); save_to_R(repmat(mouse(allSteps>0),1,2,2))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/corr_ovc_with_obj.mat', 'data')

% Fig 3G - Bayesian decoder
nbins = 100
all_error_p_env_split = nan(size(data_files,1),4);
all_error_p_env_split_rand = nan(size(data_files,1),4, nReps);
for i = 1:size(data_files,1)
      fprintf("%d/%d\n", i, size(data_files,1))
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    load(fullfile(f, 'data.mat'));
    load(fullfile(f,'celltypes.mat'));
   
    
    load(fullfile(f,'Fall.mat'));

    spks = spks(logical(iscell(:,1)),:)>0;
    
    allPCs(isnan(allPCs)) = 0;
    allOVCs(isnan(allOVCs)) = 0;
    pcs = logical(allPCs);
    ovc = logical(allOVCs);
    allPCs = ovc & ~pcs; % OVCs is any cells that are ovcs, but not pcs
    
    allLoc = objLoc'; % Use the location relative to the object
    
    mean_error_p_env = nan(4,1);
    for n = 1:2 % DOel object or not
        for s = 1:2 % MSC environmnet or not
            idx = (n-1)*2+s;
            incEnv = env==s & DO==n;
             inc = incEnv'&runFrames&incObjLoc';
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
                bin_loco = discretize(loco,nbins);
               
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
save('data/all_obj_error.mat', 'data')


all_error_p_env_rand = squeeze(nanmean(all_error_p_env_split_rand,3)) .* [2,2,2,2]; % convert to cm
datatemp = [save_to_R(all_error_p_env_rand(allSteps<0,:)); save_to_R(all_error_p_env_rand(allSteps>0,:))];
training = [save_to_R(ones(sum(allSteps<0),4)); save_to_R(ones(sum(allSteps>0),4)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4)); save_to_R(repmat(mouse(allSteps>0),1,4))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/all_obj_error_rand.mat', 'data')


%% Fig 3H - Error OVC - rand
pc_error_p_env_split = nan(size(data_files,1),4);
pc_error_p_env_rand = nan(size(data_files,1),4);


for i = 1:size(data_files,1)
      fprintf("%d/%d\n", i, size(data_files,1))
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    load(fullfile(f, 'data.mat'));
    load(fullfile(f,'celltypes.mat'));
    allPCs(isnan(allPCs)) = 0;
    
    load(fullfile(f,'Fall.mat'));
    spks = spks>0;
    
    allPCs(isnan(allPCs)) = 0;
    allOVCs(isnan(allOVCs)) = 0;
    pcs = logical(allPCs);
    ovc = logical(allOVCs);
    allPCs = ovc & ~pcs; % OVCs is any cells that are ovcs, but not pcs
    
    load(fullfile(f, 'data.mat'));
    
    allLoc = objLoc'; % Use the location relative to the object
    
    
    mean_error_p_env = nan(4,1);
    min_PCs = min(sum(allPCs));
    %% Do Bayesian decoder
    if min_PCs>0
        for n = 1:2 % DOel object or not
            for s = 1:2 % MSC environmnet or not
                idx = (n-1)*2+s;
                incEnv = env==s & DO==n;
                inc = incEnv'&runFrames&incObjLoc';
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
save('data/ovc_error.mat', 'data')