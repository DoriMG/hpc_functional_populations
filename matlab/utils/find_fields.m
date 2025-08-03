function [rmSmoothed,rmErr,logicFields,peakTrackPos,fieldWidth,fieldBins]= find_fields(ROI, binnedPos,uniqueBins, df_f, applyPCcriteria, binSize)

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

peakDeltaF=max(rmSmoothed);
peakTrackPos = find(rmSmoothed == peakDeltaF);
if length(peakTrackPos) > 0
peakTrackPos=peakTrackPos(1); %If several bins, might want to take middle one?/ Centre of mass instead of just peak value?
else
    peakTrackPos = NaN;
end

logicFields = 0;
switch applyPCcriteria
    case 1
    meanFluorPos=mean(rmSmoothed,'omitnan');
    
    findField=find(rmSmoothed > meanFluorPos); % for certain number of adjacent points. Remember indexing in 2s. Index x2 is actual track pos
    if ~isempty(findField)
         diffInd=diff(findField); %look for continuos areas of track where avg fluoresence values are above threshhold ^
         jumps=find([diffInd inf]>1); % find discontinuities/ jumps in this array(i.e. where neighbouring spatial bin is not above fluor threshold), defining the edge of the field
         widthField=diff([0 jumps]); 
         aboveThreshFieldWidth=widthField(widthField >= 5); % preliminary def of a field is 5+ continuous bins (10cm) above threshold
         if ~isempty(aboveThreshFieldWidth)
             wideEnoughFields=find(ismember(widthField,aboveThreshFieldWidth));
             logicFields=ones(1,size(wideEnoughFields,2));
            
             for fieldNum=1:size(wideEnoughFields,2)
                 wideEnoughFieldsInd=wideEnoughFields(fieldNum);
                 field=[jumps(wideEnoughFieldsInd)-(widthField(wideEnoughFieldsInd)-1):jumps(wideEnoughFieldsInd)];%indexing from original subset of values that were above thresh (not spatial bins)
                 fieldBins{fieldNum}=findField(field); %now index refers to spatial bins corresponding to each candidate field
                 fieldWidth{fieldNum}=length(fieldBins{fieldNum})*binSize;

             end
                 if ~isempty(fieldBins)
                     allFields=fieldBins{:};
                 else
                 end
         else
             logicFields=0;
             fieldWidth{1}=NaN;
             fieldBins{1}=NaN;
         end
    else
        logicFields=0;
        fieldWidth{1}=NaN;
        fieldBins{1}=NaN;
    end
      
    if sum(logicFields) >=1 % if at least one of candidate fields reach field criteria, this cell is a candidate PC
       placeCellFlag=1;     

    else
         placeCellFlag=0; 
    end

    case 0
       fieldWidth = [];
       fieldBins =[];
end
end