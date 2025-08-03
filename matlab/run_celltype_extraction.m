data_files = importdata('metadata/data_files.mat'); % list with data folders/files from the base folder
basefolder = 'D:\Frankfurt\Papers\Sussex\hippocampal_pops\v3\data';


for i = 1:size(data_files,1)
    fprintf("%d/%d\n", i, size(data_files,1))
    folder = fullfile( basefolder, strtrim(data_files{i} ));
    [f,~] = fileparts(folder);
    [allPC_ps, allSeq_ps, allPCs, allSeqs] = findCellTypes(f);
end