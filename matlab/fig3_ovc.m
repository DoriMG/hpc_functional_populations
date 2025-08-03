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
            inc = incEnv'&inc_loc&runFrames&incObjLoc';
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
diff_bins(:,:,1) = nanmean(pf_centres(:,:,1:7),3);
% familiar location control object
diff_bins(:,:,2) = nanmean(pf_centres(:,:,8:13),3);
% reward location
diff_bins(:,:,3) = nanmean(pf_centres(:,:,14:20),3);
datatemp = [save_to_R(diff_bins(allSteps<0,1:4,:)); save_to_R(diff_bins(allSteps>0,1:4,:))];
training = [save_to_R(ones(sum(allSteps<0),4,3)); save_to_R(ones(sum(allSteps>0),4,3)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,4,3)); save_to_R(repmat(mouse(allSteps>0),1,4,3))];
data = [datatemp, training(:,1), mice(:,1)];
save('data/pf_centres_binned_ovc.mat', 'data')


%% Place cell fluorescence
folders = importdata('more_folders.mat');
basefolder = 'D:\data\';
nBins = 20;
corr_p_env = nan(size(folders,1), 4);
for i = 2:size(folders,1)
    folder = fullfile( basefolder, strtrim(folders{i} ));
    [f,~] = fileparts(folder);
    
    %% Import and clean data
    allLoc = importdata(fullfile(f,'locomotion.mat'));
    Fall = importdata(fullfile(f, ['suite2p' filesep 'plane0' filesep 'Fall.mat']));
    traces = Fall.F-(0.7*Fall.Fneu);
    traces = traces(logical(Fall.iscell(:,1)),:);
    traces = findDeltaF(traces, 1:100);
    load(fullfile(f,'cellTypes_17Feb_nobehav_1ms_min.mat'))
    
    allLoc = importdata(fullfile(f,'locomotion.mat'));
    inc_loc = allLoc>50/295 & allLoc <= 250/295;
    
    velRes = getVelRegressorsBayes(allLoc, length(allLoc))*295/7.51;
    % get traversals
    tras = findTraversals(allLoc, 0.5*max(allLoc));
    
    runFrames = velRes>1;
    allLoc = allLoc*295;
    envDown = getEnvVals (f,length(allLoc));
    env = envDown(1,:);
    env = round(env);
    DO = envDown(4,:) ~= 0 & envDown(4,:) ~= 100;
    DO = DO+1; % 1 is no DOel, 2 is DOel
    load(fullfile(f,'cellTypes_17Feb_nobehav.mat'))
    allPCs(isnan(allPCs))=0;
    
    traces_norm= traces-min(traces,[],2);
    traces_norm = traces_norm./max(traces_norm, [], 2);
     load(fullfile(f,'cellTypes_23Feb_ovc.mat'));
    allPCs(isnan(allPCs)) = 0;
    ovc = logical(allPCs);
    pc_data = load(fullfile(f,'cellTypes_23Feb_pc_for_ovc_comp.mat'));
    load(fullfile(f,'cellTypes_23Feb_pc_for_ovc_comp.mat'));
    allPCs(isnan(allPCs)) = 0;
    pcs = logical(allPCs(:,1:4));
    
    locOfObj = envDown(4,:);
    locOfObj(locOfObj==0) = 100;
    objLoc = allLoc'-(locOfObj);
    % incObjLoc = abs(objLoc) <= 100;
    incObjLoc = objLoc > (-100)& objLoc < (100); % 
        inc_loc = allLoc>100 & allLoc<300; % 50 around control obj
    %% Do Bayesian decoder
    for n = 2 % DOel object or not
        for s = 1:2 % MSC environmnet or not
            idx = (n-1)*2+s;
            incEnv = env==s & DO==n;
            
            % get running periods
            ovc_cells = ovc(:,idx) & ~pcs(:,idx);
            inc = runFrames & incObjLoc & incEnv;
            spks = traces(:,inc);
            tempLoc = objLoc(inc)';
            tempLoc = tempLoc-min(tempLoc);
            tempLoc = tempLoc/max(tempLoc);
           
            [pcLocObj, pcLocSm]= placefieldsBayes(spks, round(tempLoc*nBins), nBins);
            
             inc = runFrames & inc_loc' & incEnv;
            spks = traces(:,inc);
            tempLoc = allLoc(inc)';
            tempLoc = tempLoc-min(tempLoc);
            tempLoc = tempLoc/200;
           
            [pcLoc, pcLocSm]= placefieldsBayes(spks, round(tempLoc*nBins), nBins);
            
            ovc_cells = find(ovc_cells);
            corr_temp = [];
            for c = 1:length(ovc_cells)
                incTemp = ~isnan(pcLocObj(ovc_cells(c),:)) & ~isnan(pcLoc(ovc_cells(c),:));
                if sum(incTemp)>15
                    cc = corrcoef(pcLocObj(ovc_cells(c),incTemp), pcLoc(ovc_cells(c),incTemp));
                    corr_temp = [corr_temp, cc(1,2)];
                else
                    corr_temp = NaN;
                end
            end            
            corr_p_env(i,idx) = nanmean(corr_temp);
           
        end
    end
end
% save('fluor_p_loc.mat', 'fluor_p_folder_loc', 'fluor_p_folder_loc_nopc')


%% shuffles
folders = importdata('more_folders.mat');
basefolder = 'D:\data\';
nBins = 20;
corr_p_env_rand = nan(size(folders,1), 4);
nReps = 100;
for i = 2:size(folders,1)
    folder = fullfile( basefolder, strtrim(folders{i} ));
    [f,~] = fileparts(folder);
    
    %% Import and clean data
    allLoc = importdata(fullfile(f,'locomotion.mat'));
    Fall = importdata(fullfile(f, ['suite2p' filesep 'plane0' filesep 'Fall.mat']));
    traces = Fall.F-(0.7*Fall.Fneu);
    traces = traces(logical(Fall.iscell(:,1)),:);
    traces = findDeltaF(traces, 1:100);
    load(fullfile(f,'cellTypes_17Feb_nobehav_1ms_min.mat'))
    
    allLoc = importdata(fullfile(f,'locomotion.mat'));
    inc_loc = allLoc>50/295 & allLoc <= 250/295;
    
    velRes = getVelRegressorsBayes(allLoc, length(allLoc))*295/7.51;
    % get traversals
    tras = findTraversals(allLoc, 0.5*max(allLoc));
    
    runFrames = velRes>1;
    allLoc = allLoc*295;
    envDown = getEnvVals (f,length(allLoc));
    env = envDown(1,:);
    env = round(env);
    DO = envDown(4,:) ~= 0 & envDown(4,:) ~= 100;
    DO = DO+1; % 1 is no DOel, 2 is DOel
    load(fullfile(f,'cellTypes_17Feb_nobehav.mat'))
    allPCs(isnan(allPCs))=0;
    
    traces_norm= traces-min(traces,[],2);
    traces_norm = traces_norm./max(traces_norm, [], 2);
     load(fullfile(f,'cellTypes_23Feb_ovc.mat'));
    allPCs(isnan(allPCs)) = 0;
    ovc = logical(allPCs);
    pc_data = load(fullfile(f,'cellTypes_23Feb_pc_for_ovc_comp.mat'));
    load(fullfile(f,'cellTypes_23Feb_pc_for_ovc_comp.mat'));
    allPCs(isnan(allPCs)) = 0;
    pcs = logical(allPCs(:,1:4));
    
    locOfObj = envDown(4,:);
    locOfObj(locOfObj==0) = 100;
    objLoc = allLoc'-(locOfObj);
    % incObjLoc = abs(objLoc) <= 100;
    incObjLoc = objLoc > (-100)& objLoc < (100); % 
        inc_loc = allLoc>100 & allLoc<300; % 50 around control obj

    %% Do Bayesian decoder
    for n = 2 % DOel object or not
        for s = 1:2 % MSC environmnet or not
            idx = (n-1)*2+s;
            incEnv = env==s & DO==n;
            
            % get running periods
            ovc_cells = ovc(:,idx) & ~pcs(:,idx);
            inc = runFrames & incObjLoc & incEnv;
            spks = traces(:,inc);
            tempLoc = objLoc(inc)';
            tempLoc = tempLoc-min(tempLoc);
            tempLoc = tempLoc/max(tempLoc);
           
            [pcLocObj, pcLocSm]= placefieldsBayes(spks, round(tempLoc*nBins), nBins);
            
             inc = runFrames & inc_loc' & incEnv;
            spks = traces(:,inc);
            tempLoc = allLoc(inc)';
            tempLoc = tempLoc-min(tempLoc);
            tempLoc = tempLoc/200;
           
            [pcLoc, pcLocSm]= placefieldsBayes(spks, round(tempLoc*nBins), nBins);
            
            ovc_cells = find(ovc_cells);
            corr_all_cells = nan(length(ovc_cells),1);
            for c = 1:length(ovc_cells)
                corr_temp = [];
                for rep = 1:nReps
                    rno = randi(size(pcLoc,1));
                    incTemp = ~isnan(pcLocObj(ovc_cells(c),:)) & ~isnan(pcLoc(rno,:));

                    
                if sum(incTemp)>15
                    cc = corrcoef(pcLocObj(ovc_cells(c),incTemp), pcLoc(rno,incTemp));
                    corr_temp = [corr_temp, cc(1,2)];
                else
                    corr_temp = NaN;
                end
                end
                corr_all_cells(c) = nanmean(corr_temp);
            end            
            corr_p_env_rand(i,idx) = nanmean(corr_all_cells);
           
        end
    end
end


% save('corr_p_env.mat', 'corr_p_env', 'corr_p_env_rand')
%load('corr_p_env.mat');
temp = cat(3,corr_p_env(:,3:4), corr_p_env_rand(:,3:4));
datatemp = [save_to_R(temp(allSteps<0,:,:)); save_to_R(temp(allSteps>0,:,:))];
training = [save_to_R(ones(sum(allSteps<0),2,2)); save_to_R(ones(sum(allSteps>0),2,2)*2)];
mice = [save_to_R(repmat(mouse(allSteps<0),1,2,2)); save_to_R(repmat(mouse(allSteps>0),1,2,2))];
data = [datatemp, training(:,1), mice(:,1)];
save('R/data/corr_ovc_with_obj.mat', 'data')