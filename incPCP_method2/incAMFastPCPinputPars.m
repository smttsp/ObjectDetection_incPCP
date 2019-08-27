function[pcpPars] = incAMFastPCPinputPars(dname)
%  
% Usage:
%       [pcpPars] = incAMFastPCPinputPars(name)
%
% Input:
%       name        'default',      Default parameters
%                   'none',         Create the pcpPars structure but leave it uninitialized. 
%  
%
% Output:
%       pcpPars: Structure with parameters for the incAMFastPCP method
%       
%       pcpPars.showFlag            If 1, show intermediate results.
%       pcpPars.grayFlag            If 1, treat inpt data as grayscale.
%       pcpPars.saveFlag            If 1, save low-rank and sparse estimation (as jpg files).
%       pcpPars.lambda              Regularization parameter for the incAMFastPCP method.
%       pcpPars.fullInit            If 1, use the default initialization method
%       pcpPars.stepInit            In the the default initialization method, this value represents
%                                   the (temporal) sub-sample rate between frames.
%       pcpPars.sparseGT            Sparse ground-truth
%       pcpPars.L0                  Initial low-rank approximation (currently, not in use)
%       pcpPars.S0                  Initial sparse approximation (currently, not in use)
%       pcpPars.backgroundThresh    Threshold use to determine if the background has radically changed
%       pcpPars.vminShow            Minimum value use to scale the sparse approximation only for displaying
%                                   purposes
%       pcpPars.vmaxShow            Maximum value use to scale the sparse approximation only for displaying
%                                   purposes
%       pcpPars.url                 Input data is from a streaming (live) video
%                             
%  
% Example:
%       pcpPars = incAMFastPCPInputPars('default')
%  
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


pcpPars = struct('showFlag', [], 'grayFlag', [], 'saveFlag', [], ...
                 'lambda', [], 'fullInit', [], 'stepInit', [], ...
                 'sparseGT', [], 'L0', [], 'S0', [], 'backgroundThresh',[], 'backgroundStable',[], ...
                 'url', [], 'frameNoff', [], ...
                 'TI', [], 'TImu', [], 'TIcsrLoops', [], 'TIadjustDiff', [], ...
                 'baseAlpha',[], 'baseTras', [], 'alphaThresh', [], 'cudaFlag', [], 'vecFlag', []);


                 

                 
if(nargin == 0)
  dname='none';
end


switch lower(dname)

  case{'none'}
    % do nothing

  case{'default'}

    pcpPars = setDefault_ShowSave(pcpPars);
    pcpPars = setDefault_GT(pcpPars);
    pcpPars = setDefault_fullInit_BG(pcpPars);
    
    pcpPars.lambda  = 0.025*1.5;

    pcpPars.TI      = 0;

    pcpPars.vecFlag  = 1;
    pcpPars.cudaFlag = 0;  
    
  case{'default_cuda'}

    pcpPars = setDefault_ShowSave(pcpPars);
    pcpPars = setDefault_GT(pcpPars);
    pcpPars = setDefault_fullInit_BG(pcpPars);
    
    pcpPars.lambda = 0.025*1.5;

    pcpPars.TI     =  0;

    pcpPars.vecFlag  = 0;
    pcpPars.cudaFlag = 1;  
    
    
  case{'ti_standard'}

    pcpPars = setDefault_ShowSave(pcpPars);
    pcpPars = setDefault_GT(pcpPars);
    pcpPars = setDefault_fullInit_BG(pcpPars);

    pcpPars.vecFlag  = 0;
    pcpPars.lambda   = 0.025*1.5;

    
    pcpPars.TI                  = 1;
    pcpPars.TIextraOneLoopSolve = 0;

    pcpPars = setDefault_TI_alpha(pcpPars);
    pcpPars = setDefault_TI_translation(pcpPars);

    pcpPars.TIadjustDiff      =  0.2;

    
  case{'ti_search','ti'}

    pcpPars = setDefault_ShowSave(pcpPars);
    pcpPars = setDefault_GT(pcpPars);
    pcpPars = setDefault_fullInit_BG(pcpPars);

    pcpPars.vecFlag  = 0;
    pcpPars.lambda   = 0.025*1.5;

    
    pcpPars.TI        =  1;

    pcpPars = setDefault_TI_ExtraLoopSearch(pcpPars);
    pcpPars = setDefault_TI_alpha(pcpPars);
    
    pcpPars = setDefault_TI_translation(pcpPars);

    pcpPars.TIadjustDiff      =  0.2;

    pcpPars.cudaFlag = 0;
     
  otherwise 
    disp('Unknown method.');
    pcpPars = [];


end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[pcpPars] = setDefault_ShowSave(pcpPars)

    pcpPars.grayFlag  = 0;
    pcpPars.showFlag  = 1;
    pcpPars.adaptShow = 1;
    pcpPars.url       = 0;
    pcpPars.saveFlag  = 0;
    
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[pcpPars] = setDefault_GT(pcpPars)

    pcpPars.sparseGT  = [];
    pcpPars.frameNoff =  0;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[pcpPars] = setDefault_fullInit_BG(pcpPars)

    pcpPars.fullInit = 0;
    pcpPars.stepInit = 5;
    pcpPars.backgroundThresh = 10;      % this is a real number
    pcpPars.backgroundStable = 10;      % this is number of frames

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[pcpPars] = setDefault_TI_ExtraLoopSearch(pcpPars)

  incAMFastPCPdefs;

    pcpPars.TIextraOneLoopSolve    =  1;
    pcpPars.TIextraOneLoopStrategy = INCPCP_TI_REFINEMENT_SEARCH;  % other options INCPCP_TI_REFINEMENT_STANDARD 
    pcpPars.TIAngleSearchBins      = 10;
    pcpPars.TIRotationPad          = INCPCP_TI_ROTATION_ZEROPAD;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[pcpPars] = setDefault_TI_alpha(pcpPars)

    pcpPars.baseAlpha         =  4;
    pcpPars.alphaThresh       =  0.01;
    pcpPars.TIfibonacciThresh =  0.01;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[pcpPars] = setDefault_TI_alphaCSR(pcpPars)

    pcpPars.baseAlpha         =  4;
    pcpPars.alphaThresh       =  0.05;
    pcpPars.TIfibonacciThresh =  0.05;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[pcpPars] = setDefault_TI_translation(pcpPars)

    pcpPars.baseTras          =  5;
    pcpPars.TImu              =  0.01;
    pcpPars.TIcsrThresh       =  0.75;
    pcpPars.TIcsrLoops        =  2;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[pcpPars] = setDefault_TI_ExtraLoopCSR(pcpPars)

  incAMFastPCPdefs;

    pcpPars.TIextraOneLoopSolve  =  1;
    pcpPars.TIextraOneLoopStrategy = INCPCP_TI_REFINEMENT_CSR;  % other options INCPCP_TI_REFINEMENT_STANDARD / CSR
    pcpPars.TIAngleSearchBins = 10;
    pcpPars.TIRotationPad = INCPCP_TI_ROTATION_ZEROPAD;
    
    pcpPars.TI_cbpdnascMaxIter = 20;
    pcpPars.TI_cbpdnascRho     = 5e4;
    pcpPars.TI_cbpdnascAutoRho = 1;
    pcpPars.TI_cbpdnascVerbose = 0;
    pcpPars.TI_cbpdnascLambda  = 300;
    pcpPars.TI_cbpdnascMu      = 1;

return;
