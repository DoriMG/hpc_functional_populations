function envDown = getEnvVals (f,goal)
folder = findFolders(f, '*.dat');
filename = folder{1};
% open the binary file
fid = fopen(filename);
% read all data from the file into a 5-row matrix
data = fread(fid,[8 inf],'double');
% close the file
fclose(fid);
% 
% 
% for i = 1:7
%     subplot(3,3,i); plot(envNorm(i,:));
% end

% envNorm = data - min(data,[],2);
% envNorm = envNorm./max(envNorm,[],2);
time = data(1,:)-min(data(1,:));
time = time*60*60*24;
timeRec = goal/7.51;
incTime = time<=timeRec;

envNorm = data(4:end, incTime);

runidx = 1:size(envNorm,2);                                 % Index
runidxq = linspace(min(runidx), max(runidx), goal); 
envDown = interp1(runidx, envNorm', runidxq, 'linear')';   