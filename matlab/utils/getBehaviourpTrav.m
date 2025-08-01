function [licked, nSal,nNov, displacement] =  getBehaviourpTrav(f)

files = findFolders(f,'*.dat');
filename = files{1};
rewC = [];
contC = [];
fid = fopen(filename);
data = fread(fid,[8 inf],'double');
fclose(fid);
[folder, baseFileName,~] = fileparts(filename);

time = (data(1,:) - min(data(1,:)))*(24*60*60);
loc = data(3,:);
reward = data(8,:);
novWrld = data(5,:);
novObj = data(7,:);
salWrld = data(4,:);
vel = [0,diff(loc)];
vel(vel<-10) = NaN;

lickFile =findFolders(f,'*licks.csv');
lickdata = csvread(lickFile{1});
licktime = lickdata(:,1);


%% Make sure times line up
[minValue,closestIndex] = min(abs(licktime-max(time)));
lickdata = lickdata(1:closestIndex,:);
licktime = lickdata(:,1);
licks = lickdata(:,2);

% Find the traversals
tra = 1;
tras = [tra];
for i = 1:length(loc)-1
    if loc(i+1)-loc(i) < -30
        tra = tra+1;
    end
    tras = [tras, tra];
end

thresh = nanmean(licks)*1.5;
licksbn = double(licks>thresh);
licks = interp1(licktime, licksbn', time, 'linear')';


% Known locations
rewObj = 240;
normObj = 100;


% Round locations to bin them
locInt = round(loc);
nLocs = 295;

% Variables to put results in
licked = [];
nNov = [];
nSal = [];
displacement = [];

for i = 1:max(tras)
    
    
    if (abs(nanmean(novWrld(tras==i))-1)<0.1)
        displacement = [displacement,nanmean(novObj(tras==i))-normObj];
        nNov = [nNov,1];
        for k = 1:nLocs
            curLoc = k-1;
            lickOverRew(k,i) = nanmean(licks(tras==i & locInt==curLoc));
        end
        
    else
        displacement = [displacement,0];
        nNov = [nNov,0];
        for k = 1:nLocs
            curLoc = k-1;
            lickOverRew(k,i) = nanmean(licks(tras==i & locInt==curLoc));
        end
        
    end
    if nansum(lickOverRew(rewObj-15:rewObj+5,i)>0) > 0
        licked = [licked,1];
    else
        licked = [licked,0];
    end
    if round(nanmean(salWrld(tras==i))) == 1
        nSal = [nSal,0];
    elseif round(nanmean(salWrld(tras==i))) == 2
        nSal = [nSal,1];
    else
        nSal = [nSal, NaN];
    end    
end


end