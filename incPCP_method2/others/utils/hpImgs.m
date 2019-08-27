function[Z] = hpImgs(X, flagHP);

if nargin < 2
  flagHP = 2;
end


[Nrows, Ncols, nDims] = size(X);

switch flagHP

  case{1}
  for n=1:nDims,
    Z(:,:,n) = edge(X(:,:,n),'canny');
  end
  
  case{2}
  h = [ 0 -1 0; -1 4 -1; 0 -1 0];
  for n=1:nDims,
    Z(:,:,n) = conv2(X(:,:,n), h, 'same');
    
    mask = unimodal(Z(:,:,n));
    Z(:,:,n) = mask.*Z(:,:,n);
  end
  
end
  