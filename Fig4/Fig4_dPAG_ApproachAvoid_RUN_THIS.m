%% CODE for calculating the approach/avoid states both within and between assays.
%SECTION 1 collects data from all mouse folders.
%SECTION 2 calculates the approach/avoid states within assays.
%SECTION 3 calculates the approach/avoid states between assays.
%Run sections individually (i.e., section 1 --> section 2, or section 1 --> section 3)

%% SECTION 1 (collects data from all mouse folders)

clear all

%all EPM for all dPAG mice
folders{1,1} = 'D:\dPAG_201909\EPM\230';
folders{2,1} = 'D:\dPAG_201909\EPM\816';
folders{3,1} = 'D:\dPAG_201909\EPM\825';
folders{4,1}='F:\dPAG_2019_03\3_14_2019\H15_M25_S8'; %355
folders{5,1}='F:\dPAG_2019_03\3_14_2019\H16_M59_S48'; %362
folders{6,1}='F:\dPAG_2018_08\dPAG_673\EPM\H13_M23_S23'; 
folders{7,1} = 'F:\dPAG_2018_08\dPAG_674\EPM\H13_M2_S33';

%all rat sessions for all dPAG mice
folders{1,2} = 'D:\dPAG_201909\SalineRat1\230\H11_M49_S53'; %split session - fixed-just used first session!
folders{2,2} = 'D:\dPAG_201909\SalineRat1\816';
folders{3,2} = 'D:\dPAG_201909\SalineRat1\825';
folders{4,2}='F:\dPAG_2019_03\3_20_2019\355\H10_M0_S4\'; %rat 1
folders{5,2}='F:\dPAG_2019_03\3_20_2019\362\H11_M26_S26'; %rat 1
folders{6,2}='F:\dPAG_2018_08\dPAG_673\Rat\H19_M1_S10'; %Rat 1
folders{7,2}='F:\dPAG_2018_08\dPAG_674\Rat\H18_M15_S12'; %rat 1

% First figure out 'distance' score for each mouse in EPM
for assayNum = 1
    for mouseNum = 1:size(folders,1)

        behavNum = 1;
        
        cd(folders{mouseNum,assayNum})
        load('BehaviorMS.mat','openArmIndicesMS')
        load('Tracking.mat')
        load('openTopBott.mat')
        load('PlusArmNE.mat'); load('PlusArmNW.mat'); load('PlusArmSE.mat'); load('PlusArmSW.mat');
        load('neural_data.mat')

            leftClosed = PlusArmNW(1,1);
            rightClosed = PlusArmNE(2,1);
            topOpen = PlusArmNW(1,2);
            bottomOpen = PlusArmSW(2,2);
            centerPoint = mean([topOpen,bottomOpen]);

           for idxNum = 1:length(Tracking.mouse_position) 
               if Tracking.mouse_position(idxNum,2) > PlusArmNW(2,2) & Tracking.mouse_position(idxNum,2) < PlusArmSW(1,2)
                        temp1 = (Tracking.mouse_position(idxNum,1) - leftClosed) ./ (rightClosed-leftClosed);
                        temp2 = (rightClosed - Tracking.mouse_position(idxNum,1)) ./ (rightClosed-leftClosed);
                        distFromClosed(idxNum) = min(temp1,temp2);
               else
                        temp = (abs(Tracking.mouse_position(idxNum,2) - centerPoint)) ./ (bottomOpen-topOpen);
                        distFromClosed(idxNum) = temp + .5;   
               end
           end

           for idxNum = 1:length(Tracking.mouse_positionMS) 
               if Tracking.mouse_positionMS(idxNum,2) > PlusArmNW(2,2) & Tracking.mouse_positionMS(idxNum,2) < PlusArmSW(1,2)
                        temp1 = (Tracking.mouse_positionMS(idxNum,1) - leftClosed) ./ (rightClosed-leftClosed);
                        temp2 = (rightClosed - Tracking.mouse_positionMS(idxNum,1)) ./ (rightClosed-leftClosed);
                        distFromClosedMS(idxNum) = min(temp1,temp2);
               else
                        temp = (abs(Tracking.mouse_positionMS(idxNum,2) - centerPoint)) ./ (bottomOpen-topOpen);
                        distFromClosedMS(idxNum) = temp + .5;   
               end
           end
           
           distFromClosed(find(distFromClosed < 0)) = 0;
           distFromClosed(find(distFromClosed > 1)) = 1;

           distFromClosedMS(find(distFromClosedMS < 0)) = 0;
           distFromClosedMS(find(distFromClosedMS > 1)) = 1;
           
        if openTopBott==0
            distFromClosed = 1-distFromClosed;
            distFromClosedMS = 1-distFromClosedMS;            
        end
        
        while length(distFromClosedMS) < length(neural_data.C_raw)
            distFromClosedMS = [distFromClosedMS,distFromClosedMS(end)];
        end
        while length(distFromClosedMS) > length(neural_data.C_raw)
            distFromClosedMS = distFromClosedMS(1:end-1);
        end            
        
        dataAll{behavNum,assayNum}{mouseNum} = distFromClosed;
        dataAllMS{behavNum,assayNum}{mouseNum} = distFromClosedMS;

        clearvars distFromClosed distFromClosedMS
    end
end

% First figure out 'distance' score for each mouse in Rat
for assayNum = 2%:size(folders,2)
    for mouseNum = 1:size(folders,1)

        behavNum = 1;
        
        cd(folders{mouseNum,assayNum})
        load('Tracking.mat')
        load('neural_data.mat')

        distanceMouseRat = Tracking.distanceMouseRat;
        distanceMouseRat = distanceMouseRat ./ max(distanceMouseRat); distanceMouseRat = 1-distanceMouseRat;
        distanceMouseRatMS = Tracking.distanceMouseRatMS;
        distanceMouseRatMS = distanceMouseRatMS ./ max(distanceMouseRatMS); distanceMouseRatMS = 1-distanceMouseRatMS;

        while length(distanceMouseRatMS) < length(neural_data.C_raw)
            distanceMouseRatMS = [distanceMouseRatMS;distanceMouseRatMS(end)];
        end
        while length(distanceMouseRatMS) > length(neural_data.C_raw)
            distanceMouseRatMS = distanceMouseRatMS(1:end-1);
        end            
        
        
        dataAll{behavNum,assayNum}{mouseNum} = distanceMouseRat;
        dataAllMS{behavNum,assayNum}{mouseNum} = distanceMouseRatMS;

        clearvars distFromClosed distFromClosedMS
    end
end

% Now collect the relevant behaviors for both assays

for assayNum = 1:size(folders,2)
    for mouseNum = 1:size(folders,1)
        cd(folders{mouseNum,assayNum})
        behavNum = 2;
        load('Behavior.mat','freezeIndices','escapeIndices','stretchIndices','headDipIndices','approachIndices','escapeIndices','openArmIndices')
        load('BehaviorMS.mat','freezeIndicesMS','escapeIndicesMS','stretchIndicesMS','headDipIndicesMS','approachIndicesMS','escapeIndicesMS','openArmIndicesMS')
        load('neural_data.mat')
        
        if exist('freezeIndicesMS')
        while length(freezeIndicesMS) < length(neural_data.C_raw)
            freezeIndicesMS = [freezeIndicesMS;freezeIndicesMS(end)];
        end
        while length(freezeIndicesMS) > length(neural_data.C_raw)
            freezeIndicesMS = freezeIndicesMS(1:end-1);
        end
        end

        if exist('stretchIndicesMS')
        while length(stretchIndicesMS) < length(neural_data.C_raw)
            stretchIndicesMS = [stretchIndicesMS;stretchIndicesMS(end)];
        end
        while length(stretchIndicesMS) > length(neural_data.C_raw)
            stretchIndicesMS = stretchIndicesMS(1:end-1);
        end
        end        
        
        if exist('escapeIndicesMS')
        while length(escapeIndicesMS) < length(neural_data.C_raw)
            escapeIndicesMS = [escapeIndicesMS;escapeIndicesMS(end)];
        end
        while length(escapeIndicesMS) > length(neural_data.C_raw)
            escapeIndicesMS = escapeIndicesMS(1:end-1);
        end
        end

        if exist('headDipIndicesMS')
        while length(headDipIndicesMS) < length(neural_data.C_raw)
            headDipIndicesMS = [headDipIndicesMS;headDipIndicesMS(end)];
        end
        while length(headDipIndicesMS) > length(neural_data.C_raw)
            headDipIndicesMS = headDipIndicesMS(1:end-1);
        end
        end

        if exist('openArmIndicesMS')
        while length(openArmIndicesMS) < length(neural_data.C_raw)
            openArmIndicesMS = [openArmIndicesMS;openArmIndicesMS(end)];
        end
        while length(openArmIndicesMS) > length(neural_data.C_raw)
            openArmIndicesMS = openArmIndicesMS(1:end-1);
        end
        end

        if exist('approachIndicesMS')
        while length(approachIndicesMS) < length(neural_data.C_raw)
            approachIndicesMS = [approachIndicesMS;approachIndicesMS(end)];
        end
        while length(approachIndicesMS) > length(neural_data.C_raw)
            approachIndicesMS = approachIndicesMS(1:end-1);
        end
        end
        
        
        if assayNum==1
           %behavIndicesAll = [headDipIndices';freezeIndices';openArmIndices']; 
           behavIndicesAllMS = [headDipIndicesMS';freezeIndicesMS';openArmIndicesMS']; 
           behavName = {'head dip','freeze','openArm'};
        end

        if assayNum==2
           %behavIndicesAll = [approachIndices';stretchIndices';escapeIndices';freezeIndices']; 
           behavIndicesAllMS = [approachIndicesMS';escapeIndicesMS';freezeIndicesMS';]; 
           behavName = {'approach','escape','freeze'};
        end        
        
        dataAllMS{behavNum,assayNum}{mouseNum} = behavIndicesAllMS;

        dataAllMS{behavNum+1,assayNum}{mouseNum} = behavName;

        clearvars behavIndicesAll behavIndicesAllMS freezeIndices escapeIndices stretchIndices headDipIndices approachIndices escapeIndices ...
                    openArmIndices openArmIndicesMS freezeIndicesMS escapeIndicesMS stretchIndicesMS headDipIndicesMS approachIndicesMS escapeIndicesMS
    end
end

% Now add the neural data

for assayNum = 1:size(folders,2)
    for mouseNum = 1:size(folders,1)
        cd(folders{mouseNum,assayNum})
        behavNum = 4;
        load('neural_data.mat')
        
        dataAllMS{behavNum,assayNum}{mouseNum} = neural_data.C_raw;

    end
end

% Now add the coregistration matrix
%change filepaths to each mouse's cellRegistered.mat file.

mouseNum = 1;
cd('D:\dPAG_201909\footprints\230\coreg\');
load('cellRegistered.mat')
coreg = cell_registered_struct.cell_to_index_map;
coreg = coreg(:,[2,4]);
coregAll{mouseNum} = coreg;
clearvars coreg

mouseNum = 2;
cd('D:\dPAG_201909\footprints\816\coreg');
load('cellRegistered.mat')
coreg = cell_registered_struct.cell_to_index_map;
coreg = coreg(:,[2,4]);
coregAll{mouseNum} = coreg;
clearvars coreg

mouseNum = 3;
cd('D:\dPAG_201909\footprints\825\coreg');
load('cellRegistered.mat')
coreg = cell_registered_struct.cell_to_index_map;
coreg = coreg(:,[2,4]);
coregAll{mouseNum} = coreg;
clearvars coreg

mouseNum = 4;
cd('F:\dPAG_2019_03\footprints\355\Footprints_FullOutput\coreg_addShkHab');
load('cellRegistered.mat')
coreg = cell_registered_struct.cell_to_index_map;
coreg = coreg(:,[2,3]);
coregAll{mouseNum} = coreg;
clearvars coreg

mouseNum = 5;
cd('F:\dPAG_2019_03\footprints\362\coreg_shk_hab');
load('cellRegistered.mat')
coreg = cell_registered_struct.cell_to_index_map;
coreg = coreg(:,[2,3]);
coregAll{mouseNum} = coreg;
clearvars coreg

mouseNum = 6;
cd('F:\dPAG_2018_08\dPAG_673\footprints\coreg');
load('cellRegistered.mat')
coreg = cell_registered_struct.cell_to_index_map;
coreg = coreg(:,[2,3]);
coregAll{mouseNum} = coreg;
clearvars coreg

mouseNum = 7;
cd('F:\dPAG_2018_08\dPAG_674\footprints\coreg_shk_hab');
load('cellRegistered.mat')
coreg = cell_registered_struct.cell_to_index_map;
coreg = coreg(:,[2,3]);
coregAll{mouseNum} = coreg;
clearvars coreg

% Now add the neural data -- just the coregistered cells in the correct order
    for mouseNum = 1:size(folders,1)

        tempCoreg = coregAll{mouseNum};
        
        for cellNum = 1:length(tempCoreg)
           if tempCoreg(cellNum,1)==0 | tempCoreg(cellNum,2)==0
               idxToDel(cellNum) = 1;
           else
               idxToDel(cellNum) = 0;
           end
        end
        
        tempCoreg(find(idxToDel),:) = [];
        
        dataAllMS{5,1}{mouseNum} = dataAllMS{4,1}{mouseNum}(tempCoreg(:,1),:); %for assay 1
        dataAllMS{5,2}{mouseNum} = dataAllMS{4,2}{mouseNum}(tempCoreg(:,2),:); %for assay 2

        clearvars idxToDel tempCoreg
    end
    
% Determine if mouse is moving towards or away from threat

clearvars -except folders dataAllMS dataAll coregAll

for assayNum = 1:size(folders,2)
    for mouseNum = 1:size(folders,1)
        distThreatMS = dataAllMS{1,assayNum}{mouseNum};
        distThreatMS = fillmissing(distThreatMS,'nearest');
            dataAllMS{1,assayNum}{mouseNum} = distThreatMS; %sub this in-- filled in missing values

        distThreat = dataAll{1,assayNum}{mouseNum};
        distThreat = fillmissing(distThreat,'nearest');
            dataAll{1,assayNum}{mouseNum} = distThreat; %sub this in-- filled in missing values
            
            distThreatDiffMS = diff(distThreatMS);
            distThreatDiff = diff(distThreat);
            
            for sampleNum = 4:length(distThreatDiffMS)-3
                moveThreatMS(sampleNum) = nanmean(distThreatDiffMS(sampleNum-3:sampleNum+3));
            end
            for sampleNum = 4:length(distThreatDiff)-3
                moveThreat(sampleNum) = nanmean(distThreatDiff(sampleNum-3:sampleNum+3));
            end
 
            moveThreat = [nan nan  moveThreat nan nan]; moveThreat = fillmissing(moveThreat,'nearest');
            moveThreat = zscore(fillmissing(moveThreat,'nearest'));
            moveThreatMS = [nan nan moveThreatMS nan nan]; moveThreatMS = fillmissing(moveThreatMS,'nearest');
            moveThreatMS = zscore(fillmissing(moveThreatMS,'nearest'));
            
            dirThreat = moveThreat > 0; dirThreat = double(dirThreat); temp = find(dirThreat==0); dirThreat(temp)=-1;
            dirThreatMS =  moveThreatMS > 0; dirThreatMS = double(dirThreatMS); dirThreatMS(find(dirThreatMS==0)) = -1;
            
            dataAll{6,assayNum}{mouseNum} = moveThreat;
            dataAllMS{6,assayNum}{mouseNum} = moveThreatMS;

            dataAll{7,assayNum}{mouseNum} = dirThreat;
            dataAllMS{7,assayNum}{mouseNum} = dirThreatMS;
            
            clearvars dirThreat dirThreatMS moveThreatMS moveThreat distThreatMS distThreat distThreatDiff distThreatDiffMS
    end
end

% Now calculate the approach / avoid metric

for assayNum = 1
    for mouseNum = 1:size(folders,1)
            for sampleNum = 1:length(dataAllMS{1,assayNum}{mouseNum})
                if dataAllMS{7,assayNum}{mouseNum}(sampleNum) == 1
                    temp = ((dataAllMS{1,assayNum}{mouseNum}(sampleNum)).*.9) .* dataAllMS{7,assayNum}{mouseNum}(sampleNum); %(distance * .75) .* directionMove
                elseif dataAllMS{7,assayNum}{mouseNum}(sampleNum) == -1
                    temp = (1-((dataAllMS{1,assayNum}{mouseNum}(sampleNum)).*.9)) .* dataAllMS{7,assayNum}{mouseNum}(sampleNum); %([1-distance] * .75) .* directionMove                    
                end
                
                if dataAllMS{2,assayNum}{mouseNum}(1,sampleNum) == 1
                    temp = abs(temp) .* 1.11;
                end
                if dataAllMS{2,assayNum}{mouseNum}(2,sampleNum) == 1 %if it's a freeze sample
                    temp = -1; 
                end
                               
                   appAvoid{mouseNum,assayNum}(sampleNum) = temp; clearvars temp
            end           
    end
end

for assayNum = 2%:size(folders,2)
    for mouseNum = 1:size(folders,1)
            for sampleNum = 1:length(dataAllMS{4,assayNum}{mouseNum})
                if dataAllMS{7,assayNum}{mouseNum}(sampleNum) == 1
                    temp = ((dataAllMS{1,assayNum}{mouseNum}(sampleNum)).*1) .* dataAllMS{7,assayNum}{mouseNum}(sampleNum); %(distance * .75) .* directionMove
                elseif dataAllMS{7,assayNum}{mouseNum}(sampleNum) == -1
                    temp = (1 - ((dataAllMS{1,assayNum}{mouseNum}(sampleNum)).*1)) .* dataAllMS{7,assayNum}{mouseNum}(sampleNum); %([1-distance] * .75) .* directionMove                    
                end
                
                if dataAllMS{2,assayNum}{mouseNum}(3,sampleNum) == 1 %if it's a freeze sample
                    temp = -1; 
                end
                
                    appAvoid{mouseNum,assayNum}(sampleNum) = temp; clearvars temp;               
            end

    end
end

%% SECTION 2 (calculates the approach/avoid states within assays)
% must run for separately for EPM and rat. Select which assay at top.

doPCA = 0; %PCA run for HMM, not kmeans, which uses the raw neural activity.

varThresh = 60;
bottomPC = 1;
assayToUse = 1; %'1' for EPM, '2' for Rat
numClusters = 10; %SET THE NUMBER OF CLUSTERS FOR KMEANS or STATES FOR HMM HERE

useKMeans = 1; %if '0', code runs HMM

removeNonAppAvoid = 0; %keep '0', not used.

appAvoidDel = appAvoid;

for mouseNum = 1:size(folders,1)
    for assayNum = assayToUse
            caActivity = dataAllMS{4,assayNum}{mouseNum}';

            if removeNonAppAvoid == 1
                idxDelAppAvoid{mouseNum,assayNum} = find(appAvoid{mouseNum,assayNum} == 0);
                caActivity(idxDelAppAvoid{mouseNum,assayNum},:) = [];
                appAvoidDel{mouseNum,assayNum}(idxDelAppAvoid{mouseNum,assayNum}) = [];
            end
            
            if doPCA==1
                % De-mean
                caActivity = bsxfun(@minus,caActivity,mean(caActivity));
                % Do the PCA
                [coeff,score,latent,tsquared,explained,mu] = pca(caActivity);
                %determine how many PCs to include
                temp = cumsum(explained);
                temp = find(temp > varThresh);
                numPCs = min(temp);
                caActivity = score(:,bottomPC:numPCs);
            end
            
            if useKMeans==0
            [Mu, Cov, P, Pi, LL] = hmm(caActivity, length(caActivity), numClusters, 500)
            
            for numSamples = 1:size(caActivity,1)
                for numClust = 1:size(Mu,1)
                    tempDist(numClust) = pdist2(Mu(numClust,:), caActivity(numSamples,:),'euclidean');
                end
                clusterIdx(numSamples) = find(tempDist==min(tempDist)); clearvars tempDist
            end
                statesAll{mouseNum,assayNum} = clusterIdx; clearvars clusterIdx
            end
            
            if useKMeans==1
                clusterIdx = kmeans(caActivity,numClusters,'Replicates',10,'MaxIter',1000);
                statesAll{mouseNum,assayNum} = clusterIdx; clearvars clusterIdx
            end
    end
end

% do these states differentiate the approach/avoid scores?
cntr = 1;
%labels = {'state 1','state 2'}
for assayNum = assayToUse
    for mouseNum = 1:size(folders,1)
            states = statesAll{mouseNum,assayNum};
            
            for clusterNum = 1:numClusters
                clusterIdx{clusterNum} = find(states==clusterNum);
            end
            
            for clusterNum = 1:numClusters
                appAvoid_MeanState{assayNum}(mouseNum,clusterNum) = nanmean(appAvoid{mouseNum,assayNum}(clusterIdx{clusterNum}));
            end            
    end
end
   
% do concatenated states differentiate app/avoid score?

statesConcat = [];
appAvoidConcat = [];
for assayNum = assayToUse
    for mouseNum = 1:size(folders,1)
            states = statesAll{mouseNum,assayNum};
            
            [garb, idxMax] = max(appAvoid_MeanState{assayNum}(mouseNum,:));
            [garb, idxMin] = min(appAvoid_MeanState{assayNum}(mouseNum,:));

            statesTemp = zeros(1,length(states));
            statesTemp(find(states==idxMin)) = 1;
            statesTemp(find(states==idxMax)) = 2;
            
            statesAll{mouseNum,assayNum} = statesTemp;
            
            statesConcat = [statesConcat, statesAll{mouseNum,assayNum}];
            appAvoidConcat = [appAvoidConcat, appAvoid{mouseNum,assayNum}];
    end
end

% Do stats -- is approach/avoid significantly different between states?
            State1Idx = find(statesConcat==1);
            State2Idx = find(statesConcat==2);
            appAvoid_MeanState1 = nanmean(appAvoidConcat(State1Idx));
            appAvoid_MeanState2 = nanmean(appAvoidConcat(State2Idx));
            
            [p_concat,h]=ranksum(appAvoidConcat(State1Idx),appAvoidConcat(State2Idx));

            meanConcat = [mean(appAvoidConcat(State1Idx)), mean(appAvoidConcat(State2Idx))];
            seConcat = [std(appAvoidConcat(State1Idx))./sqrt(length(appAvoidConcat(State1Idx))), std(appAvoidConcat(State2Idx))./sqrt(length(appAvoidConcat(State2Idx)))]

            figure(2002)
            bar(meanConcat); hold on;
            errorbar(meanConcat,seConcat,'LineStyle','none','Color','k')
            ylabel('avoid/approach score (-1 to 1)')
            labels = {'avoid state','approach state'};
            set(gca, 'XTickLabel', labels);
            box off;
            title(['n avoid=', num2str(length(appAvoidConcat(State1Idx))), '; n apprch=', num2str(length(appAvoidConcat(State2Idx)))])
            
% Calculate bootstrap distribution
            %bootstrap app/avoid values
            iterTotal = 100;
            
            clearvars appMeans avoidMeans
            
            for iterNum = 1:iterTotal

                     for mouseNum = 1:size(folders,1)
                        for assayNum = assayToUse
                                caActivity = dataAllMS{4,assayNum}{mouseNum}';

                                if removeNonAppAvoid == 1
                                    idxDelAppAvoid{mouseNum,assayNum} = find(appAvoid{mouseNum,assayNum} == 0);
                                    caActivity(idxDelAppAvoid{mouseNum,assayNum},:) = [];
                                    appAvoidDel{mouseNum,assayNum}(idxDelAppAvoid{mouseNum,assayNum}) = [];
                                end

                                if doPCA==1
                                    % De-mean
                                    caActivity = bsxfun(@minus,caActivity,mean(caActivity));
                                    % Do the PCA
                                    [coeff,score,latent,tsquared,explained,mu] = pca(caActivity);
                                    %determine how many PCs to include
                                    temp = cumsum(explained);
                                    temp = find(temp > varThresh);
                                    numPCs = min(temp);
                                    caActivity = score(:,bottomPC:numPCs);
                                end
                                
                                %randomize the calcium pca'ed activity
                                caActivity = caActivity(randperm(length([1:length(caActivity)])),:);                                
                                
                                if useKMeans==0
                                [Mu, Cov, P, Pi, LL] = hmm(caActivity, length(caActivity), numClusters, 400)
                                if isnan(Mu(1))
                                    A = subplot(size(folders,1),size(folders,2),idx)
                                    box off; A.XTickLabel = []; A.YTickLabel = [];
                                    if mouseNum == 1
                                        title(titleAll{assayNum})
                                    end
                                    idx = idx + 1;
                                    continue
                                end
                                clearvars clusterIdx
                                for numSamples = 1:size(caActivity,1)
                                    for numClust = 1:size(Mu,1)
                                        tempDist(numClust) = pdist2(Mu(numClust,:), caActivity(numSamples,:),'euclidean');
                                    end
                                    clusterIdx(numSamples) = find(tempDist==min(tempDist)); clearvars tempDist
                                end
                                statesAll_boot{mouseNum,assayNum} = clusterIdx; clearvars clusterIdx
                                end
                                
                                if useKMeans==1
                                    clusterIdx = kmeans(caActivity,numClusters,'MaxIter',1000);
                                    statesAll_boot{mouseNum,assayNum} = clusterIdx; clearvars clusterIdx 
                                end
                        end
                    end

                    % do these states differentiate the approach/avoid scores?
                    for assayNum = assayToUse
                        for mouseNum = 1:size(folders,1)
                                states = statesAll_boot{mouseNum,assayNum};

                                for clusterNum = 1:numClusters
                                    clusterIdx{clusterNum} = find(states==clusterNum);
                                end

                                for clusterNum = 1:numClusters
                                    appAvoid_MeanState_boot{assayNum}(mouseNum,clusterNum) = nanmean(appAvoid{mouseNum,assayNum}(clusterIdx{clusterNum}));
                                end            
                        end
                    end

                    % do concatenated states differentiate app/avoid score?
                    statesConcat_boot = [];
                    appAvoidConcat_boot = [];
                    for assayNum = assayToUse
                        for mouseNum = 1:size(folders,1)
                                states = statesAll_boot{mouseNum,assayNum};

                                [garb, idxMax] = max(appAvoid_MeanState_boot{assayNum}(mouseNum,:));
                                [garb, idxMin] = min(appAvoid_MeanState_boot{assayNum}(mouseNum,:));

                                statesTemp = zeros(1,length(states));
                                statesTemp(find(states==idxMin)) = 1;
                                statesTemp(find(states==idxMax)) = 2;

                                statesAll_boot{mouseNum,assayNum} = statesTemp;

                                statesConcat_boot = [statesConcat_boot, statesAll_boot{mouseNum,assayNum}];
                                appAvoidConcat_boot = [appAvoidConcat_boot, appAvoid{mouseNum,assayNum}];
                        end
                    end

                    % Do stats -- is approach/avoid significantly different between states?
                                State1Idx_boot = find(statesConcat_boot==1);
                                State2Idx_boot = find(statesConcat_boot==2);
                                appAvoid_MeanState1_boot = nanmean(appAvoidConcat_boot(State1Idx_boot));
                                appAvoid_MeanState2_boot = nanmean(appAvoidConcat_boot(State2Idx_boot));

                                %[p_concat,h]=ranksum(appAvoidConcat(State1Idx),appAvoidConcat(State2Idx));

                                %collect result for each iteration:
                                meanConcat_boot(iterNum,:) = [mean(appAvoidConcat_boot(State1Idx_boot)), mean(appAvoidConcat_boot(State2Idx_boot))];
                                ['Just finished iteration #' num2str(iterNum)]
            end
                    
% AND plot bootstrap distribution of mean upper and lower values with
% actual mean values (red lines)
            figure(2005)
            hist(meanConcat_boot(:,1)); hold on;
            hist(meanConcat_boot(:,2)); hold on;
            plot([meanConcat(1) meanConcat(1)],[0 30],'r'); hold on;
            plot([meanConcat(2) meanConcat(2)],[0 30],'r'); hold on;
            
            xlabel('avoid/approach score (-1 to 1), red lines actual avoid/approach vals')
            ylabel('bootstrap iter count')
            box off;
            title(['n avoid=', num2str(length(appAvoidConcat(State1Idx))), '; n apprch=', num2str(length(appAvoidConcat(State2Idx)))])

% AND AND separately plot lower and upper
            figure(2006)
            subplot(2,1,1)
            hist(meanConcat_boot(:,1)); hold on;
            plot([meanConcat(1) meanConcat(1)],[0 30],'r'); hold on;
            
            xlabel('mean avoid score (-1 to 1), red lines actual mean avoid value')
            ylabel('bootstrap iter count')
            box off;
            title(['n avoid=', num2str(length(appAvoidConcat(State1Idx)))])
            xlim([-.5 0])

            subplot(2,1,2)
            hist(meanConcat_boot(:,2)); hold on;
            plot([meanConcat(2) meanConcat(2)],[0 30],'r'); hold on;
            
            xlabel('mean approach score (-1 to 1), red lines actual mean approach value')
            ylabel('bootstrap iter count')
            box off;
            title(['n approach=', num2str(length(appAvoidConcat(State2Idx)))])
            xlim([-.2 .3])
            
%% SECTION 3 (calculates the approach/avoid states between assays)

doPCA = 0; %leave this '0'
varThresh = 60; %has no effect
bottomPC = 1; %has no effect

%1=EPM,2=Rat
trainOn = 1;
testOn = 2;

useKMeans = 1; %must be '1'.

removeNonAppAvoid = 0; %leave this '0'

numClusters = 10; %SET THE NUMBER OF CLUSTERS/STATES HERE!!!

appAvoidDel = appAvoid;

for mouseNum = 1:size(folders,1)
    for assayNum = trainOn%1:size(folders,2)
            caActivity = dataAllMS{5,assayNum}{mouseNum}';

            if removeNonAppAvoid == 1
                idxDelAppAvoid{mouseNum,assayNum} = find(appAvoid{mouseNum,assayNum} == 0);
                caActivity(idxDelAppAvoid{mouseNum,assayNum},:) = [];
                appAvoidDel{mouseNum,assayNum}(idxDelAppAvoid{mouseNum,assayNum}) = [];
            end
            
            if doPCA==1
                % De-mean
                caActivity = bsxfun(@minus,caActivity,mean(caActivity));
                % Do the PCA
                [coeff,score,latent,tsquared,explained,mu] = pca(caActivity);
                %determine how many PCs to include
                temp = cumsum(explained);
                temp = find(temp > varThresh);
                numPCs = min(temp);
                caActivity = score(:,bottomPC:numPCs);
            end
            
            if useKMeans==0
            [Mu, Cov, P, Pi, LL] = hmm(caActivity, length(caActivity), numClusters, 500)
            
            for numSamples = 1:size(caActivity,1)
                for numClust = 1:size(Mu,1)
                    tempDist(numClust) = pdist2(Mu(numClust,:), caActivity(numSamples,:),'euclidean');
                end
                clusterIdx(numSamples) = find(tempDist==min(tempDist)); clearvars tempDist
            end
                statesAll{mouseNum,assayNum} = clusterIdx; clearvars clusterIdx
            end
            
            if useKMeans==1
                [clusterIdx,C] = kmeans(caActivity,numClusters,'Replicates',10,'MaxIter',1000);
                statesAll{mouseNum,assayNum} = clusterIdx; clearvars clusterIdx
                centroidsAll{mouseNum} = C; clearvars C
            end
    end
end

%use output 'centroidsAll' variable to calculate the missing 'statesAll'
%field for the missing assay
assayNum = testOn; %# of assay number to be tested on.
for mouseNum = 1:size(folders,1)
   
    neural_temp = dataAllMS{5,assayNum}{mouseNum};

    for neuralIdx = 1:length(neural_temp)    
        for clusterNum = 1:numClusters
              tempDist(clusterNum) = pdist2(centroidsAll{mouseNum}(clusterNum,:), neural_temp(:,neuralIdx)','euclidean');
        end
        [garb,idx] = min(tempDist);
        clusterIdx(neuralIdx) = idx;
    end
        statesAll{mouseNum,assayNum} = clusterIdx'; clearvars clusterIdx
end

%now we have 'statesAll' with training and testing state id's

% do these states differentiate the approach/avoid scores?
%define approach and avoid clusters with training set here.
for assayNum = 1:size(folders,2)
    for mouseNum = 1:size(folders,1)
            states = statesAll{mouseNum,assayNum};
            
            for clusterNum = 1:numClusters
                clusterIdx{clusterNum} = find(states==clusterNum);
            end
            
            for clusterNum = 1:numClusters
                appAvoid_MeanState{assayNum}(mouseNum,clusterNum) = nanmean(appAvoid{mouseNum,assayNum}(clusterIdx{clusterNum}));
            end            
    end
end

%as a brief aside (based on convo with Jonathan), do the approach and avoid clusters have the largest and
%smallest means in the testing set?
for mouseNum = 1:size(folders,1)
    
    [garb, idxMax] =  nanmax(appAvoid_MeanState{trainOn}(mouseNum,:));
    [garb, idxMin] =  nanmin(appAvoid_MeanState{trainOn}(mouseNum,:));
    trainingMinMax(mouseNum,:) = [idxMin,idxMax];
end

%what is the actual min and max number for testing set?
for mouseNum = 1:size(folders,1)
    [garb, idxMax] =  nanmax(appAvoid_MeanState{testOn}(mouseNum,:));
    [garb, idxMin] =  nanmin(appAvoid_MeanState{testOn}(mouseNum,:));
    testingMinMax(mouseNum,:) = [idxMin,idxMax];
end


% do concatenated states differentiate app/avoid score?
statesConcat = [];
appAvoidConcat = [];
for assayNum = trainOn%1:size(folders,2)
    for mouseNum = 1:size(folders,1)
            states = statesAll{mouseNum,assayNum};
            
            [garb, idxMax] = max(appAvoid_MeanState{assayNum}(mouseNum,:));
            [garb, idxMin] = min(appAvoid_MeanState{assayNum}(mouseNum,:));

            minMaxStateIdx_Train{mouseNum} = [idxMin,idxMax];
            
            statesTemp = zeros(1,length(states));
            statesTemp(find(states==idxMin)) = 1;
            statesTemp(find(states==idxMax)) = 2;
            
            statesAll{mouseNum,assayNum} = statesTemp;
            
            statesConcat = [statesConcat, statesAll{mouseNum,assayNum}];
            appAvoidConcat = [appAvoidConcat, appAvoid{mouseNum,assayNum}];
    end
end

% Stats -- is approach/avoid significantly different between states?
            State1Idx = find(statesConcat==1);
            State2Idx = find(statesConcat==2);
            appAvoid_MeanState1 = nanmean(appAvoidConcat(State1Idx));
            appAvoid_MeanState2 = nanmean(appAvoidConcat(State2Idx));
            
            [p_concat,h]=ranksum(appAvoidConcat(State1Idx),appAvoidConcat(State2Idx));

            meanConcat = [mean(appAvoidConcat(State1Idx)), mean(appAvoidConcat(State2Idx))];
            seConcat = [std(appAvoidConcat(State1Idx))./sqrt(length(appAvoidConcat(State1Idx))), std(appAvoidConcat(State2Idx))./sqrt(length(appAvoidConcat(State2Idx)))]

            figure(2002)
            subplot(1,2,1)
            bar(meanConcat); hold on;
            errorbar(meanConcat,seConcat,'LineStyle','none','Color','k')
            %ylim([-.1 .1])
            ylabel('avoid/approach score (-1 to 1)')
            labels = {'avoid state','approach state'};
            set(gca, 'XTickLabel', labels);
            box off;
            if trainOn==2
                title(['Training Rat; n avoid=', num2str(length(appAvoidConcat(State1Idx))), '; n apprch=', num2str(length(appAvoidConcat(State2Idx)))])
            end
            if trainOn==1
                title(['Training EPM; n avoid=', num2str(length(appAvoidConcat(State1Idx))), '; n apprch=', num2str(length(appAvoidConcat(State2Idx)))])
            end

            
%calculate the above for the testing assay:

statesConcat_Test = [];
appAvoidConcat_Test = [];
for assayNum = testOn%1:size(folders,2)
    for mouseNum = 1:size(folders,1)
            states = statesAll{mouseNum,assayNum};
            
            %minMaxStateIdx_Train{mouseNum} = [idxMin,idxMax];
            
            statesTemp = zeros(1,length(states));
            statesTemp(find(states==minMaxStateIdx_Train{mouseNum}(1))) = 1;
            statesTemp(find(states==minMaxStateIdx_Train{mouseNum}(2))) = 2;
            
            statesAll{mouseNum,assayNum} = statesTemp;
            
            statesConcat_Test = [statesConcat_Test, statesAll{mouseNum,assayNum}];
            appAvoidConcat_Test = [appAvoidConcat_Test, appAvoid{mouseNum,assayNum}];
    end
end

% Stats -- is approach/avoid significantly different between states?
            State1Idx_Test = find(statesConcat_Test==1);
            State2Idx_Test = find(statesConcat_Test==2);
            appAvoid_MeanState1_Test = nanmean(appAvoidConcat_Test(State1Idx_Test));
            appAvoid_MeanState2_Test = nanmean(appAvoidConcat_Test(State2Idx_Test));
            
            [p_concat_Test,h]=ranksum(appAvoidConcat_Test(State1Idx_Test),appAvoidConcat_Test(State2Idx_Test));

            meanConcat_Test = [mean(appAvoidConcat_Test(State1Idx_Test)), mean(appAvoidConcat_Test(State2Idx_Test))];
            seConcat_Test = [std(appAvoidConcat_Test(State1Idx_Test))./sqrt(length(appAvoidConcat_Test(State1Idx_Test))), std(appAvoidConcat_Test(State2Idx_Test))./sqrt(length(appAvoidConcat_Test(State2Idx_Test)))]

            figure(2002)
            subplot(1,2,2)
            bar(meanConcat_Test); hold on;
            errorbar(meanConcat_Test,seConcat_Test,'LineStyle','none','Color','k')
            %ylim([-.1 .1])
            ylabel('avoid/approach score (-1 to 1)')
            labels = {'avoid state','approach state'};
            set(gca, 'XTickLabel', labels);
            box off;
            
            if testOn==2
                title(['Testing Rat; n avoid=', num2str(length(appAvoidConcat_Test(State1Idx_Test))), '; n apprch=', num2str(length(appAvoidConcat_Test(State2Idx_Test)))])
            end
            if testOn==1
                title(['Testing EPM; n avoid=', num2str(length(appAvoidConcat_Test(State1Idx_Test))), '; n apprch=', num2str(length(appAvoidConcat_Test(State2Idx_Test)))])
            end
            text(1.5,0,['p=', num2str(p_concat_Test)],'Color','r')

