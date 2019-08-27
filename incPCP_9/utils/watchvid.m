function[] = watchvid(V, Nrows, Ncols, nDims, jfig, myTitle, t, n)

[N frames] = size(V);

if( nargin < 8 )
  kini = 1;
  kend = frames;
else
  kini = n;
  kend = n;
end

if( nargin < 7 )
  t = -1;
if( nargin < 6 )
  myTitle = [];
  if( nargin < 5 )
    jfig = 1;
    if( nargin < 4 )
      nDims = 1;
    end
  end
end
end

if( isempty( nDims ) )
 nDims = 1;
end

vmaxShow = -1e10;
vminShow =  1e10;

for k = kini:kend

  I = V(:,k);

  vmaxShow = max([I; vmaxShow]);
  vminShow = min([I; vminShow]);


  if(nDims > 1)
    I = reshape(I, Nrows, Ncols, nDims);
  else
    I = reshape(I, Nrows, Ncols);
  end

  figure(jfig); imagesc(showNormalize(I, vminShow, vmaxShow, nDims==1)); colormap gray;
  drawnow;

  if( ~isempty(myTitle) ) title(myTitle); end

  if(t > 0) pause(t); end

end

return;

% ========================

function[y] = Normalize(x)

y = ( x - min(x(:)) ) / ( max(x(:)) - min(x(:)) );

return;

% ========================

function[y vmin vmax] = showNormalize(x, vmin, vmax, grayFlag)

if nargin < 4
  grayFlag = 0;
  if nargin < 3
    vmax = max(x(:));
    if nargin < 2
      vmin = min(x(:));
    end
  end
end

if isempty(vmax)
    vmax = max(x(:));
end

if isempty(vmin)
    vmin = min(x(:));
end

  if grayFlag
    x(1) = vmin;
    x(end) = vmax;
  end

  y = ( x - vmin ) / ( vmax - vmin );

return;
