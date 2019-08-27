clear;

filename = 'incPCP_data/sunny_youtube_stab.mp4';  
obj = VideoReader(filename);

folder = [filename 'd/'];
if exist(folder, 'dir') ~= 7
    mkdir([filename 'd/'])
    for k = 1 : obj.NumberOfFrames  %fill in the appropriate number
      this_frame = read(obj, k);
      str = sprintf('frame%05d.bmp',k);
      imwrite(this_frame, [folder str]);
    end
end

myFlags = incAMFastPCPinputPars('default'); 
myFlags.showFlag = 1; 
myFlags.grayFlag = 0; 
myFlags.saveFlag = 0; 

incrementalPCP('incPCP_data/Cloudy_MVI_40192/', 1, 2, 10, myFlags);
% 'incPCP_data/Cloudy_MVI_40192'




dctImg = [4 8 12; 16 20 24]
q = 3;
dctImgq = round(dctImg/q)*q + q/2








