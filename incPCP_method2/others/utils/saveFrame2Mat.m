function[] = saveFrame2Mat(In, k, fname, cudaFlag)

if nargin < 4
  cudaFlag = 0;
end

        if cudaFlag == 0
          I = In;
        else
          I = gather( In );
        end

        if( k<10 )
          save( sprintf('%s_frm000%d.mat', fname, k), 'I', '-v7.3' );
        end

        if( (k>=10) && (k<100) )
          save( sprintf('%s_frm00%d.mat', fname, k), 'I', '-v7.3' );
        end

        if( (k>=100) && (k<1000) )
          save( sprintf('%s_frm0%d.mat', fname, k), 'I', '-v7.3' );
        end

        if( (k>=1000) && (k<10000) )
          save( sprintf('%s_frm%d.mat', fname, k), 'I', '-v7.3' );
        end

        if( k>=10000 )
          disp('number of frames greater than 10,000 ... abort saving more frames')
          return;
        end


return;
