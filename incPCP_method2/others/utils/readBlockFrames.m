function[V] = readBlockFrames(basedir, Imgs, msize, kIni, kEnd, grayFlag, vecFlag, urlFlag, cudaFlag)


if nargin < 9
  cudaFlag = 0;
  if nargin < 8
    urlFlag = 0;
    if nargin < 7
      vecFlag = 1;
      if nargin < 6
        grayFlag = 0;
      end
    end
  end
end


  if kEnd > kIni % several images will be read
    
    if cudaFlag == 0,
      V = zeros(msize(1), kEnd-kIni+1);
    else
      V = zeros(msize(1), kEnd-kIni+1, 'gpuArray');    
    end

    for k=kIni:kEnd

      fname = get_fname(basedir, Imgs(k).name, urlFlag);
      I     = get_image(fname, grayFlag);
    
      if(vecFlag == 1) 
        if cudaFlag == 0,
          V(:, k-kIni+1) = I(:);
        else
          V(:, k-kIni+1) = gpuArray(I(:));
        end
        
      else 
        error('vecFlag must be 1 if kEnd > kIni');
      end

    end % _END_ FOR(k)

  else  % Only one image will be read
  
    k = kIni;
    
    fname = get_fname(basedir, Imgs(k).name, urlFlag);
    I     = get_image(fname, grayFlag);
  
    if cudaFlag == 1
      V = gpuArray(I);
    else
      V = vec(I, vecFlag);
    end
    
  end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[fname] = get_fname(basedir, fileName, urlFlag)

    if urlFlag == 0
      fname = strcat(basedir, fileName);
    else
      fname = basedir;
    end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[I] = get_image(fileName, grayFlag)

      if grayFlag == 1
        I = rgb2gray(double(imread(fileName))/255.0);
%          I = double(rgb2gray(imread(fileName)))/255.0;
      else
        I = double(imread(fileName))/255.0;
      end
      
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[y] = vec(x, flag)

if flag == 1
  y = x(:);
else
  y = x;
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

