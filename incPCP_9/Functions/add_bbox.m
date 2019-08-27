function [ ] = add_bbox( in, cnt)
if nargin < 2
    cnt = 1;
end

thresh = 10;
st = regionprops(in, 'BoundingBox', 'Area' );

for ii= 1 : length(st)
    Areai(ii)= st(ii).Area;
end
blobs= find(Areai==max(Areai)); 
figure(cnt), imshow(in) 
rectangle('Position',[st(blobs).BoundingBox(1),st(blobs).BoundingBox(2),...
      st(blobs).BoundingBox(3),st(blobs).BoundingBox(4)], 'EdgeColor','r','LineWidth',2 );


% st = regionprops(in, 'BoundingBox', 'Area' );
% [maxArea, indexOfMax] = max([st.Area]);
% figure(cnt), imshow(in) 
% 
% rectangle('Position',[st(indexOfMax).BoundingBox(1),st(indexOfMax).BoundingBox(2),...
%     st(indexOfMax).BoundingBox(3),st(indexOfMax).BoundingBox(4)], 'EdgeColor','r','LineWidth',2 )



end

