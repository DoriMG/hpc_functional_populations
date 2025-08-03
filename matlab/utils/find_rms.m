function [rmSmoothed,dwellMapSmoothed,rmErr]= find_rms(binnedPos,uniqueBins, df_f, binSize)

%calculate rate maps for input ROI and determine whether ROI is a place
%cell or not (optional)
%Inputs:
%ROI #
%binnedPos:
%uniqueBins: 
%df_f: 
%applyPCcriteria: 

binsInd=uniqueBins;

for i=1:size(binsInd,1)
        bin=binsInd(i);
        framePos=find(binnedPos==bin);
        dwellTime(i)=length(framePos);  
        rm(i)=mean(df_f(framePos),'omitnan');          
        rmErr(i)=std(df_f(framePos))/(sqrt(length(df_f(framePos))));
end


%%%%smooth ratemaps
kern_width = ceil(binSize * 2.5); %in bins not cm
kern=fspecial('gaussian', [1 kern_width], binSize);
rmSmoothed=imfilter(rm,kern,'same', 'conv');

dwellTime(dwellTime==0)=NaN; %set non-visited locations to NaN before smoothing. NB dwell map is not especially useful but included more to look at behaviour
dwellMapSmoothed=imfilter(dwellTime, kern, 'same', 'conv');

normMap = rmSmoothed ./ dwellMapSmoothed;
