
function[D, L, S, stats] = incrementalPCP(basedir, rank, innerLoops, winFrames, myFlags)
%  
%  [D, L] = incrementalPCP(basedir, iniK, rank, innerLoops)
%  
%  basedir      : directory where video frames are located
%  rank         : estimated rank (defaul 1)
%  innerLoops   : number of inner loops
%  winFrames    : number of frames 'to be remembered'
%  myFlags      : see incAMFastPCPinputPars.m file
%  
%  
%  
%  Examples: 
%  
%  % Standard PCP (fixed camera)
%  myFlags = incAMFastPCPinputPars('default');
%  [~, ~, ~, myStats] = incrementalPCP('./neovision3-1920x1088/', 1, 3, 50, myFlags);
%  
%  % Standard PCP (fixed camera) with CUDA-enabled functions
%  myFlags = incAMFastPCPinputPars('default'); myFlags.cudaFlag = 1;
%  [~, ~, ~, myStats] = incrementalPCP('./neovision3-1920x1088/', 1, 3, 50, myFlags);
%  
%  % Rigid transform invariant PCP (use myFlags.cudaFlag = 1 to use CUDA-enabled functions)
%  myFlags = incAMFastPCPinputPars('TI_search'); myFlags.baseTras = 10; myFlags.baseAlpha=4;
%  [~, ~, ~, myStats] = incrementalPCP('./lank3-640x480_T10-A05/', 1, 3, 50, myFlags);
%  
%  
% Authors
% =======
% 
% Paul Rodriguez   prodrig@pucp.pe
% Brendt Wohlberg  brendt@lanl.gov
% 
%  
% Related papers
% ==============
%  
%  [1] Paul Rodriguez, Brendt Wohlberg, "A Matlab Implementation of a Fast Incremental 
%      Principal Component Pursuit Algorithm for Video Background Modeling", IEEE 
%      International Conference on Image Processing (ICIP), (Paris, France), October, 2014. 
%  
%  [2] Paul Rodriguez, Brendt Wohlberg, "Incremental Principal Component Pursuit for Video 
%      Background Modeling", submitted Springer Journal of Mathematical Imaging and Vision 
%      (JMIV),  2015.
%  
%  [3] Paul Rodriguez, "Real-time Incremental Principal Component Pursuit for Video 
%      Background Modeling on the TK1", GPU Technical Conference (GTC), (San Jose, 
%      CA, USA), March, 2015.
%  
%  [4] Paul Rodriguez, Brendt Wohlberg, "Translational and Rotational Jitter Invariant 
%      Incremental Principal Component Pursuit for Video Background Modeling", accepted IEEE 
%      International Conference on Image Processing (ICIP), (Quebec, Canada), September, 2015
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
%  For other usage, please contact the authors.pcpPars
%  



% ======================
% Initial setup / checks
% ======================

addSubdirs('incrementalPCP');

grayFlag = myFlags.grayFlag;
showFlag = myFlags.showFlag;
saveFlag = myFlags.saveFlag;

lambda = myFlags.lambda;

if myFlags.cudaFlag == 1
  myFlags.vecFlag = 0;
end

vecFlag = myFlags.vecFlag;


if nargin < 5
   myFlags = incAMFastPCPinputPars('default');
   if nargin < 4
      winFrames = 30;
      if nargin < 3
        innerLoops = 1;
      end
   end
end


% Get images / video properties
[Nrows Ncols nDims frames Imgs] = getImgsProperties(basedir, grayFlag, myFlags.url);

Ndata = Nrows*Ncols*nDims;
Npix  = Nrows*Ncols;


if myFlags.url
  if frames == -1
    stats.kEnd = 10000;
  else    
    stats.kEnd = frames;
  end
else
    stats.kEnd = frames;
end


if myFlags.adaptShow
  if myFlags.cudaFlag == 0
    vmaxShow = -1e10;
    vminShow =  1e10;
  else
    vmaxShow = gpuArray(-1e10);
    vminShow = gpuArray(1e10);
  end
end

bgChangeFlag = 0;



% ===============
% Initialization
% ===============


% --- CLK is ticking (init) ---
t = tic;
% -------------------------------

if myFlags.fullInit 

  step = myFlags.stepInit;           % step (use in initialization)
  [U Sigma V] = initIncPCP(rank, step, winFrames, {Nrows, Ncols, nDims, frames, basedir, Imgs}, ...
                         {U, Sigma, V, rank}, showFlag, grayFlag, myFlags.url);

  stats.kIni = rank0+winFrames*step;

else % incremental intialization

  % read input video
  rank0 = 1;
  D = readBlockFrames(basedir, Imgs, [Nrows*Ncols*nDims, frames], 1, rank0, grayFlag, vecFlag, myFlags.url, myFlags.cudaFlag);

  
  [U Sigma] = qr(D(:), 0);  % D has rank0 columns
  
  if myFlags.cudaFlag == 0
    V = 1;
  else
    V = gpuArray(1);
  end 
  
  L = D;
  stats.kIni = rank0+1;
end

stats.Tinit = toc(t);  % initialization elapsed time


% ---------------------------



% === CLK is ticking (global) ====
t = tic;
% ================================


for k=stats.kIni:stats.kEnd, 

  t0 = tic;  % CLK per loop

  % read current frame
  curFrame = readBlockFrames(basedir, Imgs, [Nrows*Ncols*nDims, frames], k, k, grayFlag, vecFlag, myFlags.url, myFlags.cudaFlag);

  % if myFlags.TI == 0 --> alnFrame = curFrame; otherwise call TI (transform invariant)
  [alnFrame, my_h, my_hT, alphaEst] = get_TI_frame( curFrame, L, myFlags, Ncols, Nrows, nDims, Npix);
  

  %% -------------------------------------
  %% -------------------------------------

    [U Sigma V] = rank1IncSVD(U, Sigma, V, alnFrame(:), 1);  

  %% -------------------------------------
  %% -------------------------------------

  if(k > stats.kIni)
    Lold = L;
  else
    if myFlags.cudaFlag == 0
      Lold = [];
    else
      Lold = gpuArray([]);
    end
    Ldist = zeros(stats.kEnd - stats.kIni + 1, 1);
  end

  
  % ------------
  % inner loops
  % ------------
  
  for l=1:innerLoops,
  
  
    % compute current low-rank approximation ( if myFlags.TI == 0 --> L = myL = U S V')
    [L myL] = compute_LowRank(U, Sigma, V, alphaEst, my_h, Nrows, Ncols, nDims, Npix, myFlags);

    
    % compute current sparse approximation
    S = shrink( curFrame - myL, lambda);

    
    if(l==innerLoops)  break; end

    
    if(l==1)
      pFrame = alnFrame;  % 1st incSVD is applied to 'original' frame. NOTE: if myFlags.TI == 0 --> alnFrame = curFrame.
    else
      pFrame = r;
    end
    
    
    % Compute residual. if myFlags.TI == 0 --> r = curFrame-S; 
    [r, myL, my_h, my_hT, alphaEst] = computeResidual(curFrame, S, L, myL, Nrows, Ncols, nDims, Npix, myFlags, ...
                                                      my_h, my_hT, alphaEst);
    
    
    % rank-1 replace
    [U Sigma V] = rank1RepSVD(U, Sigma, V, rank, pFrame(:), r(:)); % check shkVideo idea    
    
    
  end   % _END_ innerLoops
  % ----------------------

  
  
  % ----------------------
  % -- Check BG change ---
  % ----------------------
    
  [Ldist U Sigma V bgChangeFlag vRows] = checkBGChange(k, stats.kIni, L, Lold, Ldist, U, Sigma, V, Nrows, Ncols, ...
                                                       bgChangeFlag, myFlags, curFrame);
    
  
  % ---------------
  % -- Downdate ---
  % ---------------

  if(vRows >= winFrames)    %downdate (1st col)
     [U Sigma V] = rank1DwnSVD(U, Sigma, V, 1);
  end


  % -----   CLK per loop   -----  
  stats.Tframe{k} = toc(t0);
  % ----------------------------


  % ===================================
  %    Show / Save / Compute distances
  % ===================================


  % -- Apply inverse transform (i.e. align sparse approximation for display purposes)
  Sorig = S;
  if(myFlags.TI == 1)
     S = forwardTranform( S, -alphaEst, my_hT, myFlags.alphaThresh, Nrows, Ncols, nDims, Npix, myFlags.cudaFlag);
  end
  
  
  % compute distance
  stats.dist{k} = computeDistance(S, myFlags, k, Nrows, Ncols, nDims);

  
  % compute maximum/minimum (sparse approximation, for display purposes)
  if myFlags.adaptShow
     vmaxShow = max([S(:); vmaxShow]);
     vminShow = min([S(:); vminShow]);
  end    

  
  % save data
  if saveFlag
    saveVideoFrame(S, k, './vid/S_incPCP_', Nrows, Ncols, nDims, vminShow, vmaxShow, grayFlag, myFlags.cudaFlag);
    saveVideoFrame(L, k, './vid/L_incPCP_', Nrows, Ncols, nDims, [], [], grayFlag, myFlags.cudaFlag);
    saveFrame2Mat(Sorig, k, './vid/S_incPCP_', myFlags.cudaFlag);
    
    if myFlags.url == 1
      saveVideoFrame(curFrame, k, './vid/Orig_', Nrows, Ncols, nDims, [], [], grayFlag, myFlags.cudaFlag);    
    end
  end

  
  % Show current frames
  if showFlag
    if(Ndata < 921600) % video's size less than 640x480x3
      pause(0.03);
    end
    
    
%     if k == 20
    
    S2 = reshape( S, [Nrows, Ncols, nDims]);
    S2 = S2(:,:,1)+S2(:,:,2)+S2(:,:,3);
    
    S2 = abs(S2);
    S1 = uint8(abs(S2)>0.35)*255;
   
    
%     figure(2), imshow(S2)
%     figure(3), imshow(S1)
%     
    
     S1=medfilt2(S1);
    
    sizeof = size(S2,1);
    coeff = round(sizeof/4);
    
    s_1 = S1(1:coeff,:);
    s_2 = S1(coeff+1:2*coeff,:);
    s_3 = S1(2*coeff+1:3*coeff,:);
    s_4 = S1(3*coeff+1:end,:);
    
    se1 =  strel('square',1);
    se2 =  strel('square',2);
    se3 =  strel('square',3);
    se4 =  strel('square',4);
    
    se_1 =  strel('rectangle',[2 1]);
    se_2 =  strel('rectangle',[5 2]);
    se_3 =  strel('rectangle',[8 3]);
    se_4 =  strel('rectangle',[10 4]);
    
    s_1 = imerode(s_1,se1);
    s_2 = imerode(s_2,se2);
    s_3 = imerode(s_3,se3);
    s_4 = imerode(s_4,se4);

    s_1 = imdilate(s_1,se_1);
    s_2 = imdilate(s_2,se_2);
    s_3 = imdilate(s_3,se_3);
    s_4 = imdilate(s_4,se_4);
    
    h1 = size(s_1,1);
    h2 = size(s_2,1);
    h3 = size(s_3,1);

%     close all force
%     if k == 50
        

        [group_x_left1,group_x_right1,group_y_top1,group_y_bottom1] =combine_blobs(s_1,10);
        [group_x_left2,group_x_right2,group_y_top2,group_y_bottom2] =combine_blobs(s_2,55);   
        [group_x_left3,group_x_right3,group_y_top3,group_y_bottom3] =combine_blobs(s_3,90);
        [group_x_left4,group_x_right4,group_y_top4,group_y_bottom4] =combine_blobs(s_4,150);
        

        
        x_left = [group_x_left1 group_x_left2 group_x_left3 group_x_left4];
        x_right = [group_x_right1 group_x_right2 group_x_right3 group_x_right4];
        y_top = [group_y_top1 group_y_top2 + h1 group_y_top3+ h1 + h2 group_y_top4+ h1 + h2 + h3];
        y_bottom = [group_y_bottom1 group_y_bottom2 + h1 group_y_bottom3+ h1 + h2 group_y_bottom4+ h1 + h2 + h3];
        rects = horzcat(x_left',y_top',  (x_right-x_left+1)', (y_bottom - y_top+1)');
        rects = combine_vert_box(rects);
        S1 =  vertcat(s_1,s_2,s_3,s_4);   
       
%         figure(1),imshow(S1)
        figure(1), imshow(reshape( showNormalize(curFrame), [Nrows, Ncols, nDims]));
%         hold on       
        for i =1:size(x_left,2)
            try
%         rectangle('Position',[x_left(i),y_top(i),...
%             x_right(i)-x_left(i),...
%             y_bottom(i) - y_top(i)], 'EdgeColor','r','LineWidth',2 );
            if any(rects(i))~=0
                rectangle('Position',rects(i,:), 'EdgeColor','r','LineWidth',2 );
            end
            
            catch
            end
        end
%         hold off
       
%     end    

    drawnow
%     end
%    end
    
    
    
%     for i = 1:size(group_y_top3,2)
%         for j = 1:size(group_y_top2,2)
%             if group_y_top3(i) + h1 + h2 - 1 == group_y_bottom2(j) + h1                
%                 if group_x_left3(i) <= group_x_left2(j) & group_x_right3(i) >= group_x_right2(j)
%                     disp('yes');
%                     x_left(end+1) = group_x_left3(i);
%                     x_right(end+1) = group_x_right3(i);
%                     y_top(end+1) = 
%                     y_bottom(end+1) = 
%                 end
%             end
%         end
%     end
%     
%     
% x_left = [group_x_left1 group_x_left2 group_x_left3 group_x_left4];
% x_right = [group_x_right1 group_x_right2 group_x_right3 group_x_right4];
% y_top = [group_y_top1 group_y_top2 + h1 group_y_top3+ h1 + h2 group_y_top4+ h1 + h2 + h3];
% y_bottom = [group_y_bottom1 group_y_bottom2 + h1 group_y_bottom3+ h1 + h2 group_y_bottom4+ h1 + h2 + h3];
    
    
    
    
    
    
    
    
    
    
    
   

  end


end     % _END_ FOR(k)

  
% -----   CLK (global)   -----  
stats.Tfull = toc(t);
% ----------------------------


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[localDist U Sigma V bgChangeFlag vRows] = checkBGChange(k, kIni, L, Lold, localDist, U, Sigma, V, Nrows, Ncols, bgChangeFlag, myFlags, curFrame)


  [vRows ~ ] = size(V);

  
  % -- check for changes in background
  if(k > kIni)
    if myFlags.cudaFlag == 0
      localDist(k) = sum(abs(L(:) - Lold(:))) / (Nrows*Ncols);
    else
      localDist(k) = gather(sum(abs(L(:) - Lold(:)))) / (Nrows*Ncols);      
    end
  end


  if( (k > kIni+1) )
    Lfrac = localDist(k) / localDist(k-1);
  else
    Lfrac = 1;
  end

  if( myFlags.TI == 1)
    Lfrac = Lfrac*myFlags.TIadjustDiff;  % scale down this fraction since for TI case, diferences could be due misalignement
  end

  % verbose
%    [Lfrac vRows]
  
  
  % If fraction Lfrac is greater than threshold and
  % (i) background is considered stable then re-initialize or 
  % (ii) a change has recently been detectect (this could mean that bg is changing so we 
  %      do not have to  wait until we considere it stable) then re-initialize
      
  if( (Lfrac > myFlags.backgroundThresh) & ( (vRows >= myFlags.backgroundStable) | bgChangeFlag ) )  

    [U Sigma] = qr(curFrame(:), 0);  % D has rank0 columns

    if myFlags.cudaFlag == 0
      V = 1;
    else
      V = gpuArray(1);
    end
    
    vRows = 1;
    
    bgChangeFlag = 1;

  end  
  
  if vRows >= myFlags.backgroundStable
    bgChangeFlag = 0;
  end
  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[r, myL, my_h, my_hT, alphaEst] = computeResidual(curFrame, S, L, myL, Nrows, Ncols, nDims, Npix, myFlags, my_h, my_hT, alphaEst)

    if(myFlags.TI == 1)
    
      if(myFlags.TIextraOneLoopSolve == 1)

      incAMFastPCPdefs;
      
        switch myFlags.TIextraOneLoopStrategy
        
          case{INCPCP_TI_REFINEMENT_STANDARD}
            [~, my_h, my_hT, alphaEst] = oneLoopTI_solve(curFrame-S, L, myFlags, Nrows, Ncols, nDims, Npix);
      
          case{INCPCP_TI_REFINEMENT_SEARCH}
            alphaEst = angle_refinement(curFrame, S, L, Nrows, Ncols, nDims, Npix, myFlags, my_h, alphaEst);
                      
        end % _END_ SWITCH
        
        myL = forwardTranform(L, alphaEst, my_h, myFlags.alphaThresh, Nrows, Ncols, nDims, Npix, myFlags.cudaFlag);
            
      end % _END_ IF( TIextraOneLoopSolve )
      
                                        
      r = oneLoop_IHT(curFrame, myL, S, alphaEst, myFlags.alphaThresh, my_hT, Nrows, Ncols, nDims, Npix, myFlags.cudaFlag);      
      r = L - r;
      
      
    else  %---------------------------------------

      r = curFrame-S;
            
    end  %---------------------------------------

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[alphaEst] = angle_refinement(curFrame, S, L, Nrows, Ncols, nDims, Npix, myFlags, my_h, alphaEst)

        
  pres = myFlags.TIfibonacciThresh;

  if abs(alphaEst) > 0
     alphaI = pres*floor( alphaEst / pres);
     alphaE = pres*ceil( alphaEst / pres);
  else
     alphaI = -pres/2;
     alphaE =  pres/2;
  end
        
  range = alphaE - alphaI;
  alpha = alphaI:range/(myFlags.TIAngleSearchBins-1):alphaE;
     
  if myFlags.cudaFlag == 0
    lDist = zeros(myFlags.TIAngleSearchBins,1);
  else
    lDist = zeros(myFlags.TIAngleSearchBins,1, 'gpuArray');
  end
     
     
  for k=1:myFlags.TIAngleSearchBins

      tmp   = convRotXYfull_img(L, alpha(k), 1, myFlags.cudaFlag);
      cFest = myConv3(tmp, my_h, 'same', Nrows, Ncols, nDims, Npix, myFlags.cudaFlag ) + S;

%        lDist(k) = norm( cFest(:) - curFrame(:) );
      lDist(k) = sum( ( cFest(:) - curFrame(:) ).^2 );
  end
        
  [dummy pos] = min(lDist);
  alphaEst = alpha(pos);

  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[L myL] = compute_LowRank(U, Sigma, V, alphaEst, my_h, Nrows, Ncols, nDims, Npix, myFlags) 
    
    if(myFlags.TI == 1)
      L = reshape( U*Sigma*(V(end,:)'), Nrows, Ncols, nDims);
      myL = forwardTranform(L, alphaEst, my_h, myFlags.alphaThresh, Nrows, Ncols, nDims, Npix, myFlags.cudaFlag);

    else
      L = U*Sigma*(V(end,:)');
      
      if myFlags.vecFlag == 0
         L = reshape(L, Nrows, Ncols, nDims);         
      end 
      
      myL = L;
      
    end
    
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[alnFrame, my_h, my_hT, alphaEst] = get_TI_frame( curFrame, L, myFlags, Ncols, Nrows, nDims, Npix)

  if(myFlags.TI == 1)
  
    [alnFrame, my_h, my_hT, alphaEst] = oneLoopTI_solve(curFrame, L, myFlags, Nrows, Ncols, nDims, Npix);
  
  else

    alnFrame = curFrame;

    % not used in this case
    my_h = nan;
    my_hT = nan;
    alphaEst = nan;
    
  end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[dist, recall, precision] = computeDistance(S, myFlags, k, Nrows, Ncols, nDims, verbose)

if nargin < 7
  verbose = 0;
end


  % If ground-truth given, compute distance
  if ~isempty( myFlags.sparseGT )  % sanity check

    if( ~iscell( myFlags.sparseGT ) )   % if not cell, then we assume that GT is given as MAT-files

       GT = getMatFile(myFlags.sparseGT, k);

       GT = reshape(GT, Nrows, Ncols, length(GT)/(Nrows*Ncols));

       if(myFlags.grayFlag==1) 
          GT = mean(GT,3); 
       end

       if( ~isempty(GT) )
          dist = sum( abs( S(:) - GT(:) ) ) / (Nrows*Ncols);
          if verbose==1
             [k stats.dist{k}]
          end
        else
          disp('computeDistance: GT is empty ... should not happend');
          dist = nan;
        end

        recall = nan;
        precision = nan;
        
          
    else % we have manually segmented images

      [dist recall precision ] = GTisManual(myFlags, S, Nrows, Ncols, k);

    end  % _END_ IF ~iscell( S_GT )


  else  % ground-truth was not given
  
    dist = nan;
    recall = nan;
    precision = nan;
  
  end    % _END_ IF ~isempty( S_GT )
  
  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[alnFrame, my_h, my_hT, alphaEst] = oneLoopTI_solve(shkFrame, L, myFlags, Nrows, Ncols, nDims, Npix)

    hDim = myFlags.baseTras;
    alphaBase = myFlags.baseAlpha;

    % solve current endogenous sparse representation
    [my_h my_hT] = csr4incPCP(L, shkFrame, myFlags.TImu, myFlags.TIcsrLoops, hDim, myFlags.TIcsrThresh, 0, 0, myFlags.cudaFlag);

    % --- Apply "inverse" translation ---
    tmp = myConv3( shkFrame(:), my_hT, 'same', Nrows, Ncols, nDims, Npix, myFlags.cudaFlag );


    % solve current rot via conv1D transformation (fibonacci)
    % Find alphaEst s.t. Rot(L, alphaEst) = tmp;
     myfun = @(gamma) costRotXYfull(gamma, L, tmp, 1, myFlags.cudaFlag);
     alphaEst = FSearch(myfun, -alphaBase, alphaBase, myFlags.TIfibonacciThresh);

     
    % Improve current Transformation

     if abs(alphaEst) > myFlags.alphaThresh

        % --- Apply forward rotation
        tmp = convRotXYfull_img(L, alphaEst, 1, myFlags.cudaFlag);
      
        % --- Solve (again) for translation ---
        [my_h my_hT] = csr4incPCP(tmp, shkFrame, myFlags.TImu, myFlags.TIcsrLoops, hDim, myFlags.TIcsrThresh, 0, 0, myFlags.cudaFlag);

        % --- Apply "inverse" translation ---
        tmp = myConv3( shkFrame(:), my_hT, 'same', Nrows, Ncols, nDims, Npix, myFlags.cudaFlag);


        myfun = @(gamma) costRotXYfull(gamma, L, tmp, 1, myFlags.cudaFlag);
        alphaEst = FSearch(myfun, -alphaBase, alphaBase, myFlags.TIfibonacciThresh);

        alnFrame = convRotXYfull_img(tmp, -alphaEst, 1, myFlags.cudaFlag);

     else   
      alphaEst = 0;
      alnFrame = tmp;
     end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Lshk] = forwardTranform(L, alphaEst, my_h, alphaThresh, Nrows, Ncols, nDims, Npix, cudaFlag)

    if abs(alphaEst) > alphaThresh

      % apply rigid transformation B = T(R(U))
      Lshk = convRotXYfull_img( L, alphaEst, 1, cudaFlag);    
      Lshk = myConv3(Lshk(:), my_h, 'same', Nrows, Ncols, nDims, Npix, cudaFlag);

    else
      Lshk = myConv3(L(:), my_h, 'same', Nrows, Ncols, nDims, Npix, cudaFlag);
    end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[r] = oneLoop_IHT(shkFrame, Lshk, S, alphaEst, alphaThresh, my_hT, Nrows, Ncols, nDims, Npix, cudaFlag)

      

      r = Lshk - (shkFrame - S);   % A*X - B

      if abs(alphaEst) > alphaThresh
        r = myConv3(r(:), my_hT, 'same', Nrows, Ncols, nDims, Npix, cudaFlag);     % translation
        r = convRotXYfull_img( r, -alphaEst, 1, cudaFlag);
      else
        r = myConv3(r(:), my_hT, 'same', Nrows, Ncols, nDims, Npix, cudaFlag);     % translation
      end 

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[dist] = GTisSparseApproximation(myFlags, S, Nrows, Ncols, k, oldist)


[ ismatrix( myFlags.sparseGT )  isnumeric( myFlags.sparseGT ) ]

    % ----
    if( ismatrix( myFlags.sparseGT ) & isnumeric( myFlags.sparseGT ) )
      dist = sum( abs( S - myFlags.sparseGT(:,k) ) ) / (Nrows*Ncols);

    else

      % ----
      if isdir ( myFlags.sparseGT )

      1

        GT = getMatFile(myFlags.sparseGT, k);
        if(length(GT) > 1)
          if( ~isempty(GT) )
            dist = sum( abs( S - GT ) ) / (Nrows*Ncols);
          end
        else
            dist = olddist;
        end

      end % _END_ isdir

    end % _END_ if( ismatrix( S_GT ) & isnumeric( S_GT ) )

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Fm, Recall, Precision] =  GTisManual(myFlags, S, Nrows, Ncols, k)

  myFlags.sparseGT;
  framesGT = myFlags.sparseGT{1}; 
  fnamesGT = myFlags.sparseGT{2};

  Fm = 0;
  Recall = 0;
  Precision = 0;

  if( sum( (framesGT - k - myFlags.frameNoff) == 0 ) == 1 )

    n = find( (framesGT - k - myFlags.frameNoff) == 0 );
    fname  = fnamesGT(n).name;
    
  end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[U Sigma V] = initIncPCP(k0, step, localNFrames, imgProps, cSVD, showFlag, grayFlag, urlFlag, saveFlag)

if nargin < 9
  saveFlag = 0;
  if nargin < 8
    urlFlag = 0;
    if nargin < 7
      grayFlag = 0;
    end
  end
end

Nrows   = imgProps{1};
Ncols   = imgProps{2};
nDims   = imgProps{3};
frames  = imgProps{4};
basedir = imgProps{5};
Imgs    = imgProps{6};

U     = cSVD{1};
Sigma = cSVD{2};
V     = cSVD{3};
rank  = cSVD{4};

rank0 = length(Sigma);

for k=k0:step:k0+localNFrames*step, 

  if(k>frames) break; end;

  curFrame = readBlockFrames(basedir, Imgs, [Nrows*Ncols*nDims, frames], k, k, grayFlag, urlFlag);

  if rank0 < rank
    [U Sigma V ] = rank1IncSVD(U, Sigma, V, curFrame, 0);   % rank increasing
    rank0 = rank0+1;
  else
    [U Sigma V ] = rank1IncSVD(U, Sigma, V, curFrame, 1); 
  end

  if showFlag
    L = U*Sigma*(V(end,:)');
    figure(1); imagesc( reshape( Normalize(L), [Nrows, Ncols, nDims])); colormap gray;

    if saveFlag
      saveVideoFrame(L, k, './vid/L_incPCP_', Nrows, Ncols, nDims);
    end

  end



end

% generate rank 'average' backgrounds
V = mean(V,1);
V = repmat(V, rank, 1);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[] = addSubdirs(fname)
% Add all subdirectories of the parent directory of this
% script into the path

p0 = which(fname);
K = strfind(p0, filesep);
p1 = p0(1:K(end)-1);

mypath = genpath(p1);
path(path,mypath);


%  clear p0 K mypath p1
clear p0 K mypath p1 

return;


