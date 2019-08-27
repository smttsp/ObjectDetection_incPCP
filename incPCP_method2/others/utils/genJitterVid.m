function[] = genJitterVideo(basedir, nTras, baseAlpha, grayFlag, showFlag, saveFlag, shape, firstFrame)

if nargin < 8
  firstFrame = 0;
  if nargin < 7
    shape = 'same';
    if nargin < 6
      saveFlag = 0;
      if nargin < 5
        showFlag = 1;
        if nargin < 4
          grayFlag = 0;
          if nargin < 3
            alpha = 1.5;
            if nargin < 3
              nTras = 10;
            end
          end
        end
      end
    end
  end
end

addpath utils
addpath transformations

% ---------------------------------------------------


% Get images / video properties
[Nrows Ncols nDims frames Imgs] = getImgsProperties(basedir, grayFlag);
Npix = Nrows*Ncols;


% Initial values
N = 2*nTras + 1;

% ---------------------------------------------------

switch(shape)

  case {'same'}
    Nr2 = Nrows;
    Nc2 = Ncols;

  case {'full'}
    Nr2 = Nrows+N-1;
    Nc2 = Ncols+N-1;

  otherwise
    disp('shape is not ''same'' nor ''full''... aborting');
end


% ---------------------------------------------------


for k = 1:frames,

  % random displayment
  H = zeros(N, N);

  z = abs( randn(N, N) ); 
  [~, pos] = max( z(:) );

  n = rem((pos-1), N) + 1;
  m = fix((pos-1)/N) + 1;

  H(n,m) = 1.0;

  % random rotation
  alpha = baseAlpha*2*(rand(1)-0.5);


  if( (firstFrame == 1) && (k==1) ) % Force the 1st frame to be the canonical frame

    H = zeros(N,N);
    H(nTras+1,nTras+1) = 1;

    alpha = 0;

  end


  % ---------------------------------------------------

  % Get current frame
  curFrame = genMatVid_blockFrames(basedir, Imgs, [Nrows*Ncols*nDims, frames], k, k, grayFlag);

  % Apply rigid transformation B = T(R(U))
  tmp = convRotXYfull_img( reshape(curFrame, Nrows, Ncols, nDims), alpha, 1);     % Rotation
  shkFrame = myConv3( tmp(:), H, shape, Nrows, Ncols, nDims, Npix);               % Translation

  % ---------------------------------------------------


  if showFlag == 1
    figure(1); imagesc( shkFrame ); colormap gray;
    drawnow;
  end 


  if saveFlag == 1
     saveVideoFrame(shkFrame(:), k, './vid/shk_', Nr2, Nc2, nDims);
  end


end % _END_ FOR(k)


return

