function [ out] = color_curframe( rgb, mask )
[r c]=find(mask==255);
out=rgb;
for i = 1:size(r,1)
    out(r(i),c(i),2)=255;
end
end

