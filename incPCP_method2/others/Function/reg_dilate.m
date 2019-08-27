function [ out ] = reg_dilate( im, reg)

if nargin < 2
    reg = 1;
end

h_div = floor(size(im,1)/reg);

for i = 1:reg
    out(h_div*(i-1)+1:h_div*i, :,:) = imdilate(im(h_div*(i-1)+1:h_div*i, :,:),strel('rectangle',[(i+1)*2,(i+1)*2]));
    
end
% figure,imshow(out)
end

