function incTime = getIncTime(f, inc)

folder = findFolders(f, '*.dat');
filename = folder{1};
% open the binary file
fid = fopen(filename);
% read all data from the file into a 5-row matrix
data = fread(fid,[8 inf],'double');
% close the file
fclose(fid);

time = data(1,:)-min(data(1,:));
time = time*60*60*24;

incTime = ((1:length(inc))/7.51)<=max(time);

end