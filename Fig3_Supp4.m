%Can escape or other def. behavs be predicted by rat velocity -- GLM
clear all
%all rat sessions for all dPAG mice
folders{1} = 'D:\dPAG_201909\SalineRat1\230\H11_M49_S53'; 
folders{2} = 'D:\dPAG_201909\SalineRat1\816';
folders{3} = 'D:\dPAG_201909\SalineRat1\825';
folders{4}='F:\dPAG_2019_03\3_20_2019\355\H10_M0_S4\'; 
folders{5}='F:\dPAG_2019_03\3_20_2019\358\H10_M45_S10\'; 
folders{6}='F:\dPAG_2019_03\3_20_2019\362\H11_M26_S26'; 
folders{7}='F:\dPAG_2018_08\dPAG_673\Rat\H19_M1_S10'; 
folders{8}='F:\dPAG_2018_08\dPAG_674\Rat\H18_M15_S12'; 

%% First collect the data you'll need

for mouseNum = 1:length(folders)
   
    cd(folders{mouseNum})
    load('Tracking.mat');
    load('BehaviorMS.mat');
    load('output_CNMF-E.mat','neuron');
    
    %first remove impossible rat velocities
    %Tracking.ratVelMS(find(Tracking.ratVelMS > 30)) = nan;
    %Tracking.ratVelMS = fillmissing(Tracking.ratVelMS,'nearest');
    Tracking.ratVelMS = zscore(Tracking.ratVelMS);
    
    ratVel{mouseNum} = Tracking.ratVelMS;
    escIdx{mouseNum} = escapeIndicesMS;
    strIdx{mouseNum} = stretchIndicesMS;
    frzIdx{mouseNum} = freezeIndicesMS;
    appIdx{mouseNum} = approachIndicesMS;
    
    %xdiff = jumpsmooth1D(1:length(Tracking.ratVelMS), Tracking.ratVelMS, 1, 2000, 2, 2)   
    %testThis = smoothdata(Tracking.ratVelMS,'gaussian',6); plot(Tracking.ratVelMS); hold on; plot(testThis)
    
    ratVel{mouseNum} = Tracking.ratVelMS;    
    
    while length(ratVel{mouseNum}) > length(neuron.C)
        ratVel{mouseNum} = ratVel{mouseNum}(1:end-1);
    end
    while length(escIdx{mouseNum}) > length(neuron.C)
        escIdx{mouseNum} = escIdx{mouseNum}(1:end-1);
    end
    while length(strIdx{mouseNum}) > length(neuron.C)
        strIdx{mouseNum} = strIdx{mouseNum}(1:end-1);
    end
    while length(frzIdx{mouseNum}) > length(neuron.C)
        frzIdx{mouseNum} = frzIdx{mouseNum}(1:end-1);
    end
    while length(appIdx{mouseNum}) > length(neuron.C)
        appIdx{mouseNum} = appIdx{mouseNum}(1:end-1);
    end

    while length(ratVel{mouseNum}) < length(neuron.C)
        ratVel{mouseNum} = [ratVel{mouseNum}; ratVel{mouseNum}(end)];
    end
    while length(escIdx{mouseNum}) < length(neuron.C)
        escIdx{mouseNum} = [escIdx{mouseNum}; escIdx{mouseNum}(end)];
    end
    while length(strIdx{mouseNum}) < length(neuron.C)
        strIdx{mouseNum} = [strIdx{mouseNum}; strIdx{mouseNum}(end)];
    end
    while length(frzIdx{mouseNum}) < length(neuron.C)
        frzIdx{mouseNum} = [frzIdx{mouseNum}; frzIdx{mouseNum}(end)];
    end
    while length(appIdx{mouseNum}) < length(neuron.C)
        appIdx{mouseNum} = [appIdx{mouseNum}; appIdx{mouseNum}(end)];
    end
    
end

%% concatenate across mice

frzIdxAll = [];
escIdxAll = [];
strIdxAll = [];
appIdxAll = [];
ratVelAll = [];
ratVelFrzAll = [];

for mouseNum = 1:length(folders)
    escIdxAll = [escIdxAll; escIdx{mouseNum}];
    appIdxAll = [appIdxAll; appIdx{mouseNum}];
    strIdxAll = [strIdxAll; strIdx{mouseNum}];
    ratVelAll = [ratVelAll; ratVel{mouseNum}];
end

for mouseNum = 2:length(folders)
    frzIdxAll = [frzIdxAll; frzIdx{mouseNum}];
    ratVelFrzAll = [ratVelFrzAll; ratVel{mouseNum}];  
end


%% now get coefficent for each by GLM

%escape

responseVar = escIdxAll;
predVar = ratVelAll;

            tbl = table(predVar, responseVar);
            lm = fitglm(tbl, 'responseVar~predVar');

            coeffMouse_esc = table2array(lm.Coefficients(2:end,[1,4]));



responseVar = appIdxAll;
predVar = ratVelAll;

            tbl = table(predVar, responseVar);
            lm = fitglm(tbl, 'responseVar~predVar');

            coeffMouse_app = table2array(lm.Coefficients(2:end,[1,4]));


responseVar = frzIdxAll;
predVar = ratVelFrzAll;

            tbl = table(predVar, responseVar);
            lm = fitglm(tbl, 'responseVar~predVar');

            coeffMouse_frz = table2array(lm.Coefficients(2:end,[1,4]));


%% now generate a bootstrap for each behavior and rat velocity

iterTotal = 1000;

    for iterNum = 1:iterTotal

            responseVar = escIdxAll;
            predVar = ratVelAll;

            %shuffle ratVel predictor variable
            predVar = predVar(randperm(length(predVar)));
            
            tbl = table(predVar, responseVar);
            lm = fitglm(tbl, 'responseVar~predVar');

            coeffMouse_esc_shuff(iterNum) = table2array(lm.Coefficients(2:end,1));
    end 

    
    for iterNum = 1:iterTotal

            responseVar = appIdxAll;
            predVar = ratVelAll;

            %shuffle ratVel predictor variable
            predVar = predVar(randperm(length(predVar)));
            
            tbl = table(predVar, responseVar);
            lm = fitglm(tbl, 'responseVar~predVar');

            coeffMouse_app_shuff(iterNum) = table2array(lm.Coefficients(2:end,1));
    end    

    
    for iterNum = 1:iterTotal

            responseVar = frzIdxAll;
            predVar = ratVelFrzAll;

            %shuffle ratVel predictor variable
            predVar = predVar(randperm(length(predVar)));
            
            tbl = table(predVar, responseVar);
            lm = fitglm(tbl, 'responseVar~predVar');

            coeffMouse_frz_shuff(iterNum) = table2array(lm.Coefficients(2:end,1));

    end

%%
if 0==1
   save('output_bootstrap_RatVel') 
   load('F:\dPAG_2018_08\output_bootstrap_RatVel.mat')
end

%% ...and compare bootstrap to actual vals

        temp = length(find(coeffMouse_esc(1) > coeffMouse_esc_shuff(:))) ./ iterTotal;
    
        if temp >= .95 | temp <= .05 %p-value of bootstrap and coefficient.
           mouseEsc = 1;
        else
           mouseEsc = 0;
        end
        

        temp = length(find(coeffMouse_app(1) > coeffMouse_app_shuff(:))) ./ iterTotal;
    
        if temp >= .95 | temp <= .05 %p-value of bootstrap and coefficient.
           mouseApp = 1;
        else
           mouseApp = 0;
        end


%         temp = length(find(coeffMouse_str(1) > coeffMouse_str_shuff(:))) ./ iterTotal;
%     
%         if temp >= .95 | temp <= .05 %p-value of bootstrap and coefficient.
%            mouseStr = 1;
%         else
%            mouseStr = 0;
%         end
        

        temp = length(find(coeffMouse_frz(1) > coeffMouse_frz_shuff(:))) ./ iterTotal;
    
        if temp >= .95 | temp <= .05 %p-value of bootstrap and coefficient.
           mouseFrz = 1;
        else
           mouseFrz = 0;
        end
        

%% plot value and distributions

%load('F:\dPAG_2018_08\output_Bootstrap_1000iter_ratVel_Behavs')
   
    subplot(1,3,1)    
    if coeffMouse_esc ~= 0
    %countsEsc = histcounts(coeffMouse_esc_shuff,40);
    %A = bar(countsEsc); hold on;
    %A.BarWidth = 1; box off;   
    hist(coeffMouse_esc_shuff(:)); hold on;
    plot([coeffMouse_esc(1),coeffMouse_esc(1)],[0 250],'r')
    ylim([0 300])
    title(['escape'])
    if mouseEsc == 1
        text(coeffMouse_esc(1),25,'*','Color','k','FontSize',20)
    end
    ylabel('btsrp count')
    xlabel('GLM coefficient')
    box off;
    end
    
    if coeffMouse_frz(1) ~= 0    
    subplot(1,3,2) 
    hist(coeffMouse_frz_shuff(:)); hold on;
    plot([coeffMouse_frz(1),coeffMouse_frz(1)],[0 30],'r')
    ylim([0 300])
    title(['freeze'])
    if mouseFrz == 1
        text(coeffMouse_frz(1),25,'*','Color','k','FontSize',20)
    end
    xlabel('GLM coefficient')
    end
    box off;
    
    
    if coeffMouse_app(1) ~= 0
    subplot(1,3,3) 
    hist(coeffMouse_app_shuff(:)); hold on;
    plot([coeffMouse_app(1),coeffMouse_app(1)],[0 30],'r')
    ylim([0 300])
    title(['approach'])
    if mouseApp == 1
        text(coeffMouse_app(1),25,'*','Color','k','FontSize',20)
    end
    xlabel('GLM coefficient')
    end
    box off;
    xlim([-.007 .005])
