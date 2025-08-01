function [result] = save_to_R(array)

datapoints = 1;
for i = 1:ndims(array)
    datapoints = datapoints*size(array,i);
end

result = reshape(array, [datapoints, 1]);

for i = 1:ndims(array)
    dims = size(array);
    dims(i) = 1;
    shifted = circshift(dims,-(i-1));
    els = [1:size(array,i)]';
    new_array = repmat(els, shifted);
    new_array = shiftdim(new_array,ndims(array)-(i-1));
    result = [result,reshape(new_array, [datapoints, 1])];
end