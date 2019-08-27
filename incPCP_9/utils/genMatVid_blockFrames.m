function[V] = genMatVid_blockFrames(basedir, Imgs, msize, kIni, kEnd, grayFlag, vecFlag, urlFlag)


if nargin < 8
  urlFlag = 0;
  if nargin < 7
    vecFlag = 1;
    if nargin < 6
      grayFlag = 0;
    end
  end
end

  V = zeros(msize(1), kEnd-kIni+1);

 

  for k=kIni:kEnd

    if urlFlag == 0
      fname = strcat(basedir, Imgs(k).name);
    else
      fname = basedir;
    end

    if grayFlag == 1
      I = double(rgb2gray(imread(fname)))/255;
    else
      I = double(imread(fname))/255;
    end

    if(vecFlag == 1) 
      V(:, k-kIni+1) = (I(:));
    else 
      if(kIni == kEnd) V = I;
      else V(:, k-kIni+1) = (I(:));
      end
    end

  end


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

