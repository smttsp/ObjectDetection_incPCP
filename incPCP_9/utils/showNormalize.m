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

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

