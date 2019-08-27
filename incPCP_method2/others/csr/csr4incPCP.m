function[h, hT, H, HT] = csr4incPCP(V, D, mu, loops, hDim, csrThresh, flagFreq, flagShow, cudaFlag);
%  
% Convolutional sparse representation for incPCP  
%  
%  solves 0.5 || h*u - b ||_2 + \lambda || h ||_1
%  
% Authors
% =======
% 
% Paul Rodriguez   prodrig@pucp.pe
% Brendt Wohlberg  brendt@lanl.gov
% 
%  
% Legal
% =====
% 
%  There is currently a patent pending that covers the incremental PCP method and
%  applications that is implemented in this source code.
%  
%  For non-commercial or academic use the source code of this program can be distributed 
%  and/or modified under the terms of the GNU Lesser General Public License (LGPL) version 3 
%  as published by the Free Software Foundation (http://opensource.org/licenses/lgpl-3.0.html).
%  
%  For other usage, please contact the authors.
%  

if nargin < 9
  cudaFlag = 0;
  if nargin < 8
    flagShow = 0;
    if nargin < 7
      flagFreq = 0;
      if nargin < 6
        csrThresh = 0.75;
      end
    end
  end
end

if isempty(csrThresh)
  csrThresh = 0.75;
end

[Nrows Ncols nDims] = size(D);

if nDims == 1

  Vf = fft2(V - mean(V(:)));
  Df = fft2(D - mean(D(:)));

else 

  tmp = rgb2gray(V);
  Vf = fft2(tmp - mean(tmp(:)));

  tmp = rgb2gray(D);
  Df = fft2(tmp - mean(tmp(:)));

end


Vfc = conj(Vf);

Num = Vfc.*Df;
Dem = Vfc.*Vf;

if cudaFlag == 0
  localMu = mu*max(abs(Df(:)));
else
  localMu = gpuArray(mu)*max(abs(Df(:)));
end

% --- init ---

G = Num./(Dem + localMu);
h = real(ifft2(G, 'symmetric'));
[mask, thresh] = unimodal(h, csrThresh);  % 'adaptive shrinkage' NOTE: mask is GPU if h is.
h = h.*mask;


for k = 2:loops,

  if( sum( abs(h(:)) > 0 ) == 1), break; end

  H = fft2(h);

  G = (Num + localMu*H)./(Dem + localMu);
  h = real(ifft2(G, 'symmetric'));
  [mask, thresh] = unimodal(h, csrThresh);  % 'adaptive shrinkage'
  h = h.*mask;


end
%  

H = fft2(h) / sum(h(:));    % needed?

tmp = fftshift(h);


if rem(Nrows,2) == 0
  Nr2 = ceil( Nrows/2 );
else
  Nr2 = floor( Nrows/2 );
end

if rem(Ncols,2) == 0
  Nc2 = ceil( Ncols/2 );
else
  Nc2 = floor( Ncols/2 );
end



h = tmp(Nr2-hDim+1:Nr2+hDim+1, Nc2-hDim+1:Nc2+hDim+1);
k=1;
while( sum(h(:)) == 0 )
  h = tmp(Nr2-hDim+1-k:Nr2+hDim+1+k, Nc2-hDim+1-k:Nc2+hDim+1+k);
  k = k+1;
end

if( sum( abs(h(:)) > 0 ) > 1 )

%    disp('cond');
  [v1 p1] = max(h);
  [v2 p2] = max(v1);
  
  h(p1(p2), p2) = v2;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = h/sum(h(:));
hT = rot90(rot90(h));

if flagFreq == 1

  if cudaFlag == 0
    tmp = zeros(Nrows, Ncols);
  else
    tmp = zeros(Nrows, Ncols, 'gpuArray');
  end
  
  [hrows hcols] = size(h);
  tmp(1:hrows, 1:hcols) = h;
  H = fft2(tmp);


  if cudaFlag == 0
    tmp = zeros(Nrows, Ncols);
  else
    tmp = zeros(Nrows, Ncols, 'gpuArray');
  end

  [hrows hcols] = size(hT);
  tmp(1:hrows, 1:hcols) = hT;
  HT = fft2(tmp);

else

  H  = 0;
  HT = 0;
end


%%%%%%%%%%%%%%%%%
if flagShow
  
  if nDims == 1
    Vt = conv2(V, h, 'same');
  else

    Vt(:,:,1) = conv2(V(:,:,1), h, 'same');
    Vt(:,:,2) = conv2(V(:,:,2), h, 'same');
    Vt(:,:,3) = conv2(V(:,:,3), h, 'same');

  end


  figure(1); imagesc( Vt );  colormap gray; axis image
  figure(2); imagesc( D );  colormap gray; axis image
  figure(3); imagesc( D-Vt );  colormap gray; axis image
  drawnow;
end
