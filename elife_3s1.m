%are the coregistered cells actually the same cells -- this tries to
%determine this using their PNR and mean peak size across coregistered
%sessions.

clear all; close all;

filepath{8,2} = 'F:\dPAG_2019_03\footprints\362\';
filepath{7,2} = 'F:\dPAG_2019_03\footprints\358\Footprints_FullOutput\'; %shock
filepath{6,2} = 'F:\dPAG_2019_03\footprints\355\Footprints_FullOutput\';
filepath{4,2} = 'F:\dPAG_2018_08\dPAG_673\footprints\'; 
filepath{5,2} = 'F:\dPAG_2018_08\dPAG_674\footprints\'; 

filepath{1,2} = 'D:\dPAG_201909\footprints\230'; 
filepath{2,2} = 'D:\dPAG_201909\footprints\816'; 
filepath{3,2} = 'D:\dPAG_201909\footprints\825'; 

numNeighbors = 10;

for contextIdx = 2  
    for mouseIdx = 1:size(filepath,1)
    
            cd(filepath{mouseIdx, contextIdx})
            load('PNRmetrics.mat'); load('PF_Info.mat'); 
            
            for sessionIdx = 1:size(CellReg_ShockControl,2)
                for cellIdx = 1:length(CellReg_ShockControl(:,1))
                    if CellReg_ShockControl(cellIdx,sessionIdx) > 0
                        PNR_AllSess_shock{mouseIdx}(cellIdx,sessionIdx) = PNR_neuronAll{sessionIdx}(CellReg_ShockControl(cellIdx,sessionIdx));
                    else
                        PNR_AllSess_shock{mouseIdx}(cellIdx,sessionIdx) = nan;
                    end
                    
                    if CellReg_ShockControl(cellIdx,sessionIdx) > 0
                        peakMean_AllSess_shock{mouseIdx}(cellIdx,sessionIdx) = peakMean_neuronAll{sessionIdx}(CellReg_ShockControl(cellIdx,sessionIdx));
                    else
                        peakMean_AllSess_shock{mouseIdx}(cellIdx,sessionIdx) = nan;
                    end

                end
                    %if CellReg_ShockControl(cellIdx,sessionIdx) > 0

                        for refIdx=1:size(center_neuronAll{sessionIdx},1) %the reference neuron
                            for compIdx=1:size(center_neuronAll{sessionIdx},1) %the comparison neuron (all are compared)
                                dist(refIdx, compIdx) = pdist([center_neuronAll{sessionIdx}(refIdx,:); center_neuronAll{sessionIdx}(compIdx,:)],'euclidean');
                            end
                        end

                        %sort indices and distances by ascending order
                        for i = 1:size(dist,1)
                            [b, sortOutput] = sort(dist(i,:), 'ascend');
                            distIndicesSorted(i,:) = sortOutput;
                            distScoresSorted(i,:) = b;
                        end

                            distIndicesSorted(:,1) = [];
                            distIndicesSorted = distIndicesSorted(:,1:numNeighbors);  %keep the closest neighbors

                            distIndices_AllSess_shock{sessionIdx, mouseIdx} = distIndicesSorted;

                            clearvars distIndicesSorted distScoresSorted sortOutput b dist
                    %end

            end
    end
end

%put the distance indices into CellReg_ShockControl format
for contextIdx = 2   
    for mouseIdx = 1:size(filepath,1)
        
            cd(filepath{mouseIdx, contextIdx})
            load('PNRmetrics.mat'); load('PF_Info.mat');

            for sessionIdx = 1:size(CellReg_ShockControl,2)
                    for cellIdx = 1:length(CellReg_ShockControl(:,1))

                        if CellReg_ShockControl(cellIdx,sessionIdx) > 0
                            [none,tempIdx] = ismember(distIndices_AllSess_shock{sessionIdx,mouseIdx}(CellReg_ShockControl(cellIdx,sessionIdx),:), CellReg_ShockControl(:,sessionIdx));   
                            Dist_AllSess_Shock{mouseIdx}{cellIdx,sessionIdx} = tempIdx;
                        else
                            Dist_AllSess_Shock{mouseIdx}{cellIdx,sessionIdx} = nan;
                        end
                    end
            end
    end
end

% Use Dist_AllSess_Shock to build cell matrix of all possible values for
% PNR and meanPeak values, to be used to shuffle, for the closest neighboring cells.

for mouseIdx = 1:size(Dist_AllSess_Shock,2)
    
            cd(filepath{mouseIdx, contextIdx})
            load('PNRmetrics.mat'); load('PF_Info.mat');

   for sessionIdx = 1:size(Dist_AllSess_Shock{mouseIdx},2)
       for cellIdx = 1:size(Dist_AllSess_Shock{mouseIdx},1)
           if ~isnan(Dist_AllSess_Shock{mouseIdx}{cellIdx,sessionIdx})
            idxDel = find(Dist_AllSess_Shock{mouseIdx}{cellIdx,sessionIdx} == 0); Dist_AllSess_Shock{mouseIdx}{cellIdx,sessionIdx}(idxDel) = [];
            %PNR_AllSess_shuffle_shock{mouseIdx}{cellIdx,sessionIdx} = PNR_AllSess_shock{mouseIdx}(CellReg_ShockControl(Dist_AllSess_Shock{mouseIdx}{cellIdx,sessionIdx},sessionIdx));

            PNR_AllSess_shuffle_shock{mouseIdx}{cellIdx,sessionIdx} = PNR_AllSess_shock{mouseIdx}(Dist_AllSess_Shock{mouseIdx}{cellIdx,sessionIdx},sessionIdx);
            
            %peakMean_AllSess_shuffle_shock{mouseIdx}{cellIdx,sessionIdx} = peakMean_AllSess_shock{mouseIdx}(CellReg_ShockControl(Dist_AllSess_Shock{mouseIdx}{cellIdx,sessionIdx},sessionIdx));
            peakMean_AllSess_shuffle_shock{mouseIdx}{cellIdx,sessionIdx} = peakMean_AllSess_shock{mouseIdx}(Dist_AllSess_Shock{mouseIdx}{cellIdx,sessionIdx},sessionIdx);

           else
            PNR_AllSess_shuffle_shock{mouseIdx}{cellIdx,sessionIdx} = nan;
            peakMean_AllSess_shuffle_shock{mouseIdx}{cellIdx,sessionIdx} = nan; 
           end
       end
   end
end

%% combine across mice

%take only rat1 and shock sessions
firstSess = [4 4 4 3 3 3 3 3]; %all rat 1 sessions
secondSess = [2 2 2 2 2 2 2 2]; %all EPM sessions

for mouseNum = 1:length(filepath)
   peakMean_AllSess_shock{mouseNum} = peakMean_AllSess_shock{mouseNum}(:,[firstSess(mouseNum),secondSess(mouseNum)]);
   PNR_AllSess_shock{mouseNum} = PNR_AllSess_shock{mouseNum}(:,[firstSess(mouseNum),secondSess(mouseNum)]);
   peakMean_AllSess_shuffle_shock{mouseNum} = peakMean_AllSess_shuffle_shock{mouseNum}(:,[firstSess(mouseNum),secondSess(mouseNum)]);
   PNR_AllSess_shuffle_shock{mouseNum} = PNR_AllSess_shuffle_shock{mouseNum}(:,[firstSess(mouseNum),secondSess(mouseNum)]);
end


peakMean_AllSess_shockALL = zeros(1, size(peakMean_AllSess_shock{1},2));
for i = 1:size(peakMean_AllSess_shock,2)
        peakMean_AllSess_shockALL = [peakMean_AllSess_shockALL; peakMean_AllSess_shock{1,i}];
end
peakMean_AllSess_shockALL(1,:) = [];

% peakMean_AllSess_controlALL = zeros(1, size(peakMean_AllSess_control{1},2));
% for i = 1:size(peakMean_AllSess_control,2)
%         peakMean_AllSess_controlALL = [peakMean_AllSess_controlALL; peakMean_AllSess_control{1,i}];
% end
% peakMean_AllSess_controlALL(1,:) = [];

PNR_AllSess_shockALL = zeros(1, size(PNR_AllSess_shock{1},2));
for i = 1:size(PNR_AllSess_shock,2)
        PNR_AllSess_shockALL = [PNR_AllSess_shockALL; PNR_AllSess_shock{1,i}];
end
PNR_AllSess_shockALL(1,:) = [];

% PNR_AllSess_controlALL = zeros(1, size(PNR_AllSess_control{1},2));
% for i = 1:size(PNR_AllSess_control,2)
%         PNR_AllSess_controlALL = [PNR_AllSess_controlALL; PNR_AllSess_control{1,i}];
% end
% PNR_AllSess_controlALL(1,:) = [];

PNR_AllSess_shuffle_shockALL = cell(1, size(PNR_AllSess_shuffle_shock{1},2));
for i = 1:size(PNR_AllSess_shuffle_shock,2)
        PNR_AllSess_shuffle_shockALL = [PNR_AllSess_shuffle_shockALL; PNR_AllSess_shuffle_shock{1,i}];
end
PNR_AllSess_shuffle_shockALL(1,:) = [];

peakMean_AllSess_shuffle_shockALL = cell(1, size(peakMean_AllSess_shuffle_shock{1},2));
for i = 1:size(peakMean_AllSess_shuffle_shock,2)
        peakMean_AllSess_shuffle_shockALL = [peakMean_AllSess_shuffle_shockALL; peakMean_AllSess_shuffle_shock{1,i}];
end
peakMean_AllSess_shuffle_shockALL(1,:) = [];


%% SAME AS THE SECTION ABOVE, BUT CONSTRAIN SHUFFLE TO CLOSEST NEURONS

%bc dealing only with columns 4 and 5, delete rows that have NaN in col 4
%or 5 from all analyzed matrices

temp = isnan(PNR_AllSess_shockALL(:,1)) | isnan(PNR_AllSess_shockALL(:,2)); idxDel = find(temp == 1);

PNR_AllSess_shockALL(idxDel,:) = [];
peakMean_AllSess_shockALL(idxDel,:) = [];
PNR_AllSess_shuffle_shockALL(idxDel,:) = [];
peakMean_AllSess_shuffle_shockALL(idxDel,:) = [];

PNR_r_val = corr(PNR_AllSess_shockALL(:,1), PNR_AllSess_shockALL(:,2),'Type','Spearman','Rows','complete');
peakMean_r_val = corr(peakMean_AllSess_shockALL(:,1), peakMean_AllSess_shockALL(:,2),'Type','Spearman','Rows','complete');

%Create a bootstrap distribution to compare with these r-values.

samples = 1000;

for i = 1:samples
    for sessionIdx = 1:size(PNR_AllSess_shockALL,2)
        for cellIdx = 1:size(PNR_AllSess_shockALL,1)
               if ~isnan(PNR_AllSess_shuffle_shockALL{cellIdx,sessionIdx})
                   randOrder = randperm(numNeighbors-3); 
                   if length(PNR_AllSess_shuffle_shockALL{cellIdx,sessionIdx}) == numNeighbors
                       PNR_AllSess_shockALL_shuffleNeighbor(cellIdx,sessionIdx) = PNR_AllSess_shuffle_shockALL{cellIdx,sessionIdx}(randOrder(1));
                       peakMean_AllSess_shockALL_shuffleNeighbor(cellIdx,sessionIdx) = peakMean_AllSess_shuffle_shockALL{cellIdx,sessionIdx}(randOrder(1));
                   else
                       PNR_AllSess_shockALL_shuffleNeighbor(cellIdx,sessionIdx) = PNR_AllSess_shuffle_shockALL{cellIdx,sessionIdx}(1);
                       %PNR_AllSess_shockALL_shuffleNeighbor(cellIdx,sessionIdx) = nan;
                       peakMean_AllSess_shockALL_shuffleNeighbor(cellIdx,sessionIdx) = peakMean_AllSess_shuffle_shockALL{cellIdx,sessionIdx}(1);
                       %peakMean_AllSess_shockALL_shuffleNeighbor(cellIdx,sessionIdx) = nan;
                   end
               else
                   PNR_AllSess_shockALL_shuffleNeighbor(cellIdx,sessionIdx) = nan;
                   peakMean_AllSess_shockALL_shuffleNeighbor(cellIdx,sessionIdx) = nan;                   
               end
        end
    end
                    PNR_r_val_bootstr_neighbor(i) = corr(PNR_AllSess_shockALL(:,1), PNR_AllSess_shockALL_shuffleNeighbor(:,2),'Type','Spearman','Rows','complete');
                    peakMean_r_val_bootstr_neighbor(i) = corr(peakMean_AllSess_shockALL(:,1), peakMean_AllSess_shockALL_shuffleNeighbor(:,2),'Type','Spearman','Rows','complete');               

    clearvars PNR_AllSess_shockALL_shuffleNeighbor peakMean_AllSess_shockALL_shuffleNeighbor
end

temp = find(PNR_r_val_bootstr_neighbor > PNR_r_val);
pval_PNR_neighbor = length(temp)/samples;

temp = find(peakMean_r_val_bootstr_neighbor > peakMean_r_val);
pval_peakMean_neighbor = length(temp)/samples;

%%
figure(1)
subplot(1,2,1)
hist(PNR_r_val_bootstr_neighbor,20); hold on;
%ylim([0 1800]); xlim([.1 .7]);
arrow([PNR_r_val,150],[PNR_r_val,100],'Width', .00001, 'Length', 20, 'Color', [0 .3 .7])
title('Peak-to-noise ratio bootstrap distribution, 1000 samples')
xlabel('correlation of PNR (Spearman), rat 1 and EPM');
ylabel('count');
box off;

text(.15, 125, ['true r=' num2str(round(PNR_r_val,3))])
text(.15, 105, ['p=', num2str(round(pval_PNR_neighbor,3))])


subplot(1,2,2)
hist(peakMean_r_val_bootstr_neighbor,20); hold on;
%ylim([0 1800]); xlim([.1 .7]);
arrow([peakMean_r_val,150],[peakMean_r_val,100],'Width', .00001, 'Length', 20, 'Color', [0 .3 .7])
title('Peak mean bootstrap distribution, 1000 samples')
xlabel('correlation of mean peak (Spearman), rat 1 and EPM');
ylabel('count');
box off;

text(.15, 125, ['true r=' num2str(round(peakMean_r_val,3))])
text(.15, 105, ['p=', num2str(round(pval_peakMean_neighbor,3))])



