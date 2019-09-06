myDir = 'C:\Users\17029\Desktop\gfortran\datfiles'; %gets directory
myFiles = dir(fullfile(myDir,'*.dat')); %gets all txt files in struct


for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
  fullFileName = fullfile(myDir, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  num = importdata(fullFileName);   %or readtable
  
  dimension = length(num(1,:)) - 1;
  
  figure(k)
  
  for i = 2: (dimension + 1)
  plot(num(:,1), num(:,i))
  hold on 
  end
  
  title(num(1,[2:length(num(1,:))]))
  xlabel("Mu value")
  ylabel("Lambda Values")
  
 
end