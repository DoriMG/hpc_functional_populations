
%% function written by Kira, August 2017
%function to normalise data by a specified baseline period - 
%finds delta F/F

%send in the data to norm, and the frames within which to normalise it 
%sends out the normalised data 

function [dataNorm]=findDeltaF(data2norm, frames2norm)

%determine if data is a cell (i.e. when averaging across ALL sessions)
if iscell(data2norm)==0;
    
    %find out if the data has 2 or 3 dimensions
    t=ndims(data2norm);
    
    %for 2D data, just need to loop through 1st dim, and norm
    if t==2
        for i=1:size(data2norm,1);
            %subtract baseline, then divide by it... 
            dataNorm(i,:)= (data2norm(i,:)-(nanmean(data2norm(i,frames2norm(1):frames2norm(2)))))...
                ./nanmean(data2norm(i,frames2norm(1):frames2norm(2)));
        end
        
        %for 3D data, need to loop through 1st 2 dims and norm
    elseif t==3
        %normalise data by dividing by the mean of the normalisation period
        for i=1:size(data2norm,1);
            for j=1:size(data2norm,2);
                dataNorm(i,j,:)= (data2norm(i,j,:)-(nanmean(squeeze(data2norm(i,j,frames2norm(1):frames2norm(2))))))...
                    ./nanmean(squeeze(data2norm(i,j,frames2norm(1):frames2norm(2))));
            end
        end
    end
    
    %averaging across ALL sessions
elseif iscell(data2norm)==1;
    
    for a=1:size(data2norm,2)
        
        data_ttt=data2norm{1,a};
        
        if size(isnan(data_ttt),1)>1
            
            %find out if the data has 2 or 3 dimensions
            t=ndims(data_ttt);
            
            %for 2D data, just need to loop through 1st dim, and norm
            if t==2
                for i=1:size(data2norm,1);
                    dataNorm{1,a}(i,:)= (data_ttt(i,:)-(nanmean(squeeze(data_ttt(i,frames2norm(1):frames2norm(2))))))...
                        ./nanmean(squeeze(data_ttt(i,frames2norm(1):frames2norm(2))));
                end
                
                %for 3D data, need to loop through 1st 2 dims and norm
            elseif t==3
                %normalise data by dividing by the mean of the normalisation period
                for i=1:size(data_ttt,1);
                    for j=1:size(data_ttt,2);
                         dataNorm{1,a}(i,j,:)= (data_ttt(i,j,:)-(nanmean(squeeze(data_ttt(i,j,frames2norm(1):frames2norm(2))))))...
                            ./nanmean(squeeze(data_ttt(i,j,frames2norm(1):frames2norm(2))));
                    end
                end
            end
                      
        else
            
            dataNorm{1,a}=NaN;
            
        end
        
    end
    
end

end