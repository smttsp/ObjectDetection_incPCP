function[] = addSubdirsPath(fname)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

