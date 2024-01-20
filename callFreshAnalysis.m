%this code will parse threough an entire pre-processed data set and call freshAnalysis to analyze each cell event pulse

% ch_height = 30; De_np = 15.97;wC = 8;thresholds = [135, 1100]; filepath ='/Users/taylorthomsen/Library/CloudStorage/GoogleDrive-taylorthomsen@berkeley.edu/My Drive/Sohn lab/stahl/081623FIRSTTRYONMZ/rawDataDone/adipocytes-stahl-081623-70strainer-8V-P75-8um-try3-dev10a/0417preproc.mat'

function [finalOut] = ...
    callFreshAnalysis(filepath, ch_height, De_np, wC, thresholds)

%% 
    if nargin < 5 || isempty(thresholds)
        thresholds = [135, 1100];
        fprintf('Auto thresholds set to %3.2e, %3.2e\n',thresholds); 
    end

%%
load(filepath, 'Rf','yas2det','sampleRate')

window=movmean(yas2det, 10000); % heavy movmean to get general peak of an entire event
[~,locs]=findpeaks(window,'minPeakHeight',10000,'minPeakDistance',14000);

finalOut = ones(length(locs),20);
cellCount = 0;

figure
axes
for i = 1:length(locs)

   cellCount = cellCount + 1;
  

   pk = locs(i);
   startidx = pk-7000;
   
   endidx = pk+10000;


   
 

   int = yas2det(startidx:endidx);


   filteredData=movmean(int, 80);
    
   diff1 = diff(filteredData); %difference signal of filtered and indexed data
   
  % need to add discard or keep event
   thresholds = thresholdingTT(diff1,thresholds, cellCount, startidx);
   
   userInputAnalyze = input("Do you want to analyze? (1 for Y/2 for N)");
   if userInputAnalyze == 2
       continue;
   end

   [OUT_array] = ...
        freshAnalysis(Rf,yas2det, startidx, endidx, sampleRate, ch_height, De_np, wC, thresholds, false, false);
   
   finalOut(i,:) = OUT_array; %need to keep old rows and just append new one for each iteration

end


rows_to_remove = all(finalOut == 1, 2);
finalOut(rows_to_remove, :) = [];
% wCCol=ones(size(finalOut))*wC;
% finalOut = finalOut + wCCol;

