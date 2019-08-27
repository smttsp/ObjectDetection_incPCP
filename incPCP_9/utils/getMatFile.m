function[GT] = getMatFile(basedir, k)

  Imgs   = dir(strcat(basedir,'*.mat'));

  if length(Imgs) == 0
    disp('NOTE:');
    disp(sprintf('  There are no MAT file in directory %s (aborting) ...',basedir));

    GT  = [];
    return;
  end

  fname = strcat(basedir, Imgs(k).name);
  
  try
    GT = struct2array( load(fname) );
  catch
    GT = -1;
  end

return;


