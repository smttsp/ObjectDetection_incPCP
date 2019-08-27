function [ out ] = my_kmeans( im, k)

if nargin < 2
    k=2;
end

[r,c,d]=size(im);
out=kmeans(double(im(:)), k, 'start', 'uniform', 'emptyaction','singleton');

one=(sum(out==1));
two=(sum(out==2));

if one > two %class 1 is = 0, class 2 = 1
    out = 255*uint8(out-1);
else         %class 1 is = 1, class 2 = 0
    out = uint8(abs(out-2))*255;
end
out = reshape(out,[r,c,d]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   kmeans image segmentation
%
%   Input:
%          ima: grey color image
%          k: Number of classes
%   Output:
%          mu: vector of class means 
%          mask: clasification image mask
%
%   Author: Jose Vicente Manjon Herrera
%    Email: jmanjon@fis.upv.es
%     Date: 27-08-2005
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if nargin < 2
%     k=2;
% end
% 
% % check image
% im=double(im);
% copy=im;         % make a copy
% im=im(:);       % vectorize ima
% mi=min(im);      % deal with negative 
% im=im-mi+1;     % and zero values
% 
% s=length(im);
% 
% % create image histogram
% 
% m=max(im)+1;
% h=zeros(1,m);
% hc=zeros(1,m);
% 
% for i=1:s
%   if(im(i)>0) h(im(i))=h(im(i))+1;end;
% end
% ind=find(h);
% hl=length(ind);
% 
% % initiate centroids
% 
% out=(1:k)*m/(k+1);
% 
% % start process
% 
% while(true)
%   
%   oldmu=out;
%   % current classification  
%  
%   for i=1:hl
%       c=abs(ind(i)-out);
%       cc=find(c==min(c));
%       hc(ind(i))=cc(1);
%   end
%   
%   %recalculation of means  
%   
%   for i=1:k, 
%       a=find(hc==i);
%       out(i)=sum(a.*h(a))/sum(h(a));
%   end
%   
%   if(out==oldmu) break;end;
%   
% end
% 
% % calculate mask
% s=size(copy);
% mask=zeros(s);
% for i=1:s(1),
% for j=1:s(2),
%   c=abs(copy(i,j)-out);
%   a=find(c==min(c));  
%   mask(i,j)=a(1);
% end
% end
% 
% out=out+mi-1;   % recover real range
% 



end

