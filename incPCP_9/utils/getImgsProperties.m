function[Nrows Ncols nDims frames Imgs] = getImgsProperties(basedir, grayFlag, urlFlag);

if nargin < 3
  urlFlag = 0;
  if nargin < 2
    grayFlag = 0;
  end
end

% --------------------
% >>> Read dataset <<<

if urlFlag == 0

  Imgs   = dir(strcat(basedir,'*.jpg'));
  if length(Imgs) == 0
    Imgs   = dir(strcat(basedir,'*.bmp'));
    if length(Imgs) == 0
      Imgs   = dir(strcat(basedir,'*.png'));
    end
  end

  if length(Imgs) == 0
    disp('NOTE:');
    disp(sprintf('  There are no jpg/bmp/png images in directory %s (aborting) ...',basedir));

    Nrows = -1; Ncols = -1; nDims = -1; frames = -1;
    Imgs  = '';
    return;
  end


  frames = length(Imgs);
  fname  = strcat(basedir, Imgs(1).name);

else    % ----

  fname = basedir;  % assumed to be an url address
  frames = -1;
  Imgs = [];

end


I = imread(fname) ;
if grayFlag == 1
  I = rgb2gray(I);
end
[Nrows Ncols nDims] = size(I);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


