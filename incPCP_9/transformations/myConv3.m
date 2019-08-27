function[X] = myConv3(In, H, shape, Nrows, Ncols, nDims, Npix, cudaFlag)
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

if nargin < 8
  cudaFlag = 0;
end

  switch shape

    case{'full'}
      [Lrows, Lcols] = size(H);
      if cudaFlag == 0
        X = zeros(Nrows+Lrows-1, Ncols+Lcols-1, nDims);
      else
        X = zeros(Nrows+Lrows-1, Ncols+Lcols-1, nDims, 'gpuArray');
      end

    case{'same'}
      if cudaFlag == 0
        X = zeros(Nrows, Ncols, nDims);
      else
        X = zeros(Nrows, Ncols, nDims, 'gpuArray');
      end

    otherwise
      X = 0;    % this will force an error

  end

  
  for l=1:nDims,
    X(:,:,l) = conv2( reshape( In(1+(l-1)*Npix:Npix+(l-1)*Npix), Nrows, Ncols), H, shape);
  end


return
