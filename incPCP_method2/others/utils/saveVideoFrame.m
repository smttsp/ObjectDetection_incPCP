function[] = saveVideoFrame(In, k, fname, Nrows, Ncols, nDims, vmin, vmax, grayFlag, cudaFlag)

if nargin < 10
  cudaFlag = 0;
  if nargin < 9
    grayFlag = 0;
    if nargin < 8
      vmax = max(In(:));
      if nargin < 7
        vmin = min(In(:));
      end
    end
  end
end

if isempty(vmax)
    vmax = max(In(:));
end

if isempty(vmin)
    vmin = min(In(:));
end


        if cudaFlag == 0
          I = showNormalize(reshape( In, [Nrows, Ncols, nDims]), vmin, vmax, grayFlag);
        else
          I = gather( showNormalize(reshape( In, [Nrows, Ncols, nDims]), vmin, vmax, grayFlag) );
        end
        
        if( k<10 )
          imwrite( I, sprintf('%s_img000%d.jpg', fname, k), 'jpeg', 'Quality', 85);
        end

        if( (k>=10) && (k<100) )
          imwrite( I, sprintf('%s_img00%d.jpg', fname, k), 'jpeg', 'Quality', 85);
        end

        if( (k>=100) && (k<1000) )
          imwrite( I, sprintf('%s_img0%d.jpg', fname, k), 'jpeg', 'Quality', 85);
        end

        if( (k>=1000) && (k<10000) )
          imwrite( I, sprintf('%s_img%d.jpg', fname, k), 'jpeg', 'Quality', 85);
        end

        if( k>=10000 )
          disp('number of frames greater than 10,000 ... abort saving more frames')
          return;
        end


return;
