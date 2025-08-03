function all_mi = findMutualInfo(traces, pos_dat, track_length, bin_sz)

[~,binnedPos,binsUnique]=histnd2(pos_dat, track_length, bin_sz);
all_mi = [];
nQuart = 4;


for c = 1:size(traces,1)    
    cell_temp = traces(c,:);
    [cell,E] = discretize(cell_temp,nQuart);
    
    pi = nan(100,1);
    pj = nan(nQuart,1);
    pij = nan(100,nQuart);
    for i = 1:length(binsUnique)
        pi(i) = sum(binnedPos == i)/length(binnedPos);
        for j = 1:nQuart
            pj(j) = sum(cell == j)/length(cell);
            pij(i,j) = sum(binnedPos == i & cell' == j)/length(cell);
        end
    end
    
    mi = 0;
    for i = 1:length(binsUnique)
        for j = 1:nQuart
            mi_temp = pij(i,j)*log2(pij(i,j)/(pi(i)*pj(j)));
            if ~isnan(mi_temp)
                mi = mi+mi_temp;
            end
        end
    end
    all_mi = [all_mi, mi];
end
end