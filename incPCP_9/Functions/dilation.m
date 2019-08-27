function [ out ] = dilation( in,se )
if nargin < 2
    se = strel('rectangle',[3,3]);
end
out = imdilate(in,se);

end

