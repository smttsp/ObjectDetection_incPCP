
function[V] = genMatVid(basedir, Imgs, msize, grayFlag)
%  
%  Generate a matrix that represents a video. Each column of the matrix represents
%  a frame. The frames (images) are located in the directory 'basedir'
%  

  if nargin < 4
    grayFlag = 0;
  end

  V = zeros(msize(1), msize(2));

  for k=1:msize(2)

    fname = strcat(basedir, Imgs(k).name);
    if grayFlag == 1
      I = double(rgb2gray(imread(fname)))/255;
    else
      I = double(imread(fname))/255;
    end

    V(:, k) = (I(:));

  end

return;

