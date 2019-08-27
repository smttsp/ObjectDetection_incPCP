function[cost] = costRotXYfull(alpha, I, Ir, flagCenter, cudaFlag)
%  
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

if nargin < 5
  cudaFlag = 0;
end

[Nrows Ncols nDims] = size(I);

if nDims == 3
  Ig  = mean(I,3);
  Igr = mean(Ir,3);
%    Ig  = rgb2gray(I);
%    Igr = rgb2gray(Ir);
else
  Ig  = I;
  Igr = Ir;
end

vmin = min(Ig(:));
vmax = max(Ig(:));


if cudaFlag == 1
  Ix = zeros([Nrows Ncols nDims], 'gpuArray');
  Iy = zeros([Nrows Ncols nDims], 'gpuArray');
  Ir = zeros([Nrows Ncols nDims], 'gpuArray');
else
  Ix = zeros([Nrows Ncols nDims]);
  Iy = zeros([Nrows Ncols nDims]);
  Ir = zeros([Nrows Ncols nDims]);
end

a = pi*alpha/180;


%%%%%%%%%%%% 1st pass %%%%%%%%%%%
% X-dir

if flagCenter == 1
  Nr2 = Nrows/2;
  y = -Nr2+1:Nr2;
  
else
  y = 0:Nrows-1;
end

x = 0:Ncols/2;

if cudaFlag == 0
  wx = -1i*2*pi*x/Ncols;
  delta = -y(:)*tan(a/2);
else
  wx    = gpuArray(-1i*2*pi*x/Ncols);
  delta = gpuArray(-y(:)*tan(a/2));
end


W = exp( delta*wx );

Wx = [W, conj(W(:, end-1:-1:2)) ] ;

% apply in X

If = fft(Ig, [], 2);
Ifx = If.*Wx;
Ix =  ifft( Ifx, [], 2, 'symmetric');

%%%%%%%%%%%% 2nd pass %%%%%%%%%%%

% Y-dir
y = 0:Nrows/2;

if flagCenter == 1
  Nc2 = Ncols/2;
  x = -Nc2+1:Nc2;
else
  x = 0:Ncols-1;
end



if cudaFlag == 0
  wy = -1i*2*pi*y/Ncols;
  delta = x*sin(a);
else
  wy    = gpuArray(-1i*2*pi*y/Ncols);
  delta = gpuArray( x*sin(a) );
end

W = exp( wy(:)*delta );
Wy = [W; conj(W(end-1:-1:2,:)) ];


% apply in Y

If  = fft(Ix, [], 1);
Ify = If.*Wy;
Iy = ifft( Ify, [], 1, 'symmetric' );


%%%%%%%%%%%% 3rd pass %%%%%%%%%%%
% X-dir

% apply in X

If = fft(Iy, [], 2);
Ifx = If.*Wx;
Iest =  ifft( Ifx, [], 2, 'symmetric');

m1 = Iest < vmin;
m2 = Iest > vmax;

% adjust dynamic range
Iest = Iest.*(1-m1).*(1-m2) + vmin*m1 + vmax*m2;

% compute cost
cost = norm(Iest(:) - Igr(:)) / norm(Igr(:));

