function[] = watchvid2(V1, V2, Nrows, Ncols, nDims, jfig, myTitle, t, n)

[N frames] = size(V1);

if( nargin < 9 )
  kini = 1;
  kend = frames;
else
  kini = n;
  kend = n;
end

if( nargin < 8 )
  t = -1;
if( nargin < 7 )
  myTitle = [];
  if( nargin < 6 )
    jfig = 1;
    if( nargin < 3 )
      nDims = 1;
    end
  end
end
end

if( isempty( nDims ) )
 nDims = 1;
end


for k = kini:kend

  I1 = V1(:,k);
  I1 = reshape(I1, Nrows, Ncols, nDims);

  I2 = V2(:,k);
  I2 = reshape(I2, Nrows, Ncols, nDims);

  if(nDims == 3)

    I2(1,1,1) = 0;  I2(1,1,2) = 0; I2(1,1,3) = 0;
    I2(end,end,1) = 1;  I2(end,end,2) = 1; I2(end,end,3) = 1;

  elseif(nDims == 1)

    I2(1,1) = 0;
    I2(end, end) = 1;
  end


  I = [I1 I2];

  figure(jfig); imagesc(Normalize(I)); colormap gray; drawnow;

  if( ~isempty(myTitle) ) title(myTitle); end

  if(t > 0) pause(t); end

end



% ========================

function[y] = Normalize(x)

y = ( x - min(x(:)) ) / ( max(x(:)) - min(x(:)) );


