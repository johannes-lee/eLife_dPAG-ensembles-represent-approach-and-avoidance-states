%% Calculate many variables within Rat and EPM assays to get correlations
% SMP, updated 07/23/20

clear all

sess{1} = ''; %\Rat1'


cd('');
load('');

celltypes_all = {celltype355, celltype362, celltype816, ...
    celltype825, celltype230, celltype358, celltype674};
keepinds_all = {keepind355, keepind362, keepind816, ...
    keepind825, keepind230, keepind358, keepind674};

% Storage arrays
% EPM
epm_ca_count_ = zeros(size(sess));
epm_oa_count_ = zeros(size(sess));
epm_ca_count_p_ = zeros(size(sess));
epm_oa_count_p_ = zeros(size(sess));
epm_zc_neie_mean_ = zeros(size(sess));
epm_zc_opne_mean_ = zeros(size(sess));
epm_zc_csde_mean_ = zeros(size(sess));
epm_oa_zc_neie_mean_ = zeros(size(sess));
epm_oa_zc_opne_mean_ = zeros(size(sess));
epm_oa_zc_csde_mean_ = zeros(size(sess));
epm_ca_zc_neie_mean_ = zeros(size(sess));
epm_ca_zc_opne_mean_ = zeros(size(sess));
epm_ca_zc_csde_mean_ = zeros(size(sess));
epm_frz_percent_ = zeros(size(sess));
epm_str_percent_ = zeros(size(sess));
epm_frz_count_ = zeros(size(sess));
epm_str_count_ = zeros(size(sess));
epm_vel_mean_ = zeros(size(sess));


% RAT
rat_celltypes_all = cell(length(sess), 1);
nea_all = cell(length(sess), 1);
far_all = cell(length(sess), 1);
rat_nea_percent_ = zeros(length(sess), 1);
rat_far_percent_ = zeros(length(sess), 1);
rat_nea_count_ = zeros(length(sess), 1);
rat_far_count_ = zeros(length(sess), 1);
% for each ensemble
% neither
rat_zc_neie_mean_ = zeros(size(sess));
% open
rat_zc_opne_mean_ = zeros(size(sess));
% closed
rat_zc_csde_mean_ = zeros(size(sess));
% near for each ensemble
rat_nea_zc_neie_mean_ = zeros(size(sess));
rat_nea_zc_opne_mean_ = zeros(size(sess));
rat_nea_zc_csde_mean_ = zeros(size(sess));
% far for each ensemble
rat_far_zc_neie_mean_ = zeros(size(sess));
rat_far_zc_opne_mean_ = zeros(size(sess));
rat_far_zc_csde_mean_ = zeros(size(sess));
% Percent Behaviors
rat_frz_percent_ = zeros(size(sess));
rat_esc_percent_ = zeros(size(sess));
rat_str_percent_ = zeros(size(sess));
rat_apr_percent_ = zeros(size(sess));
% Count Behaviors
rat_frz_count_ = zeros(size(sess));
rat_esc_count_ = zeros(size(sess));
rat_str_count_ = zeros(size(sess));
rat_apr_count_ = zeros(size(sess));
% Time in each EPM region
epm_t_o_p_ = zeros(size(sess));
epm_t_c_p_ = zeros(size(sess));
epm_t_m_p_ = zeros(size(sess));

t_cells = zeros(size(sess));
n_opne = zeros(size(sess));
n_csde = zeros(size(sess));


for jj = 1:length(sess)
    cd('D:\dPAG');
    cd(sess{jj});
    load('cellRegistered');
    
    animal = sess{jj};
    
    celltypes = celltypes_all{jj};
    keepind = transpose(keepinds_all{jj});
    
    t_cells(jj) = sum(keepind);
    n_opne(jj) = sum(celltypes == 'o'&keepind == 1);
    n_csde(jj) = sum(celltypes == 'c'&keepind == 1);
    
    if (jj < 2 | jj > 5) % mice 674, 355, 358, 362
        cor_rat = cell_registered_struct.cell_to_index_map(:, 3);
    else % mice 230, 816, 825
        cor_rat = cell_registered_struct.cell_to_index_map(:, 4);
    end
    cor_epm = cell_registered_struct.cell_to_index_map(:, 2);
    cor_ref = cell_registered_struct.cell_to_index_map(:, 1);
    
    cd('EPM');
    load('Tracking.mat');
    try
        load('BehaviorMS.mat');
    catch
        load('output_CNMF-E.mat');
    end
    
    load('neural_data.mat');
    
    % Organize the data a bit:
    % Time in seconds
    timestamp = readtable('timestamp.dat');
    t_sec = table2array(timestamp(:, 'sysClock'));
    t_sec = (t_sec/1000);
    if(t_sec(size(timestamp, 1), 1) == 0)
        timestamp([size(timestamp, 1)], :) = [];
        t_sec = table2array(timestamp(:, 'sysClock'));
        t_sec = (t_sec/1000);
    end
    if(t_sec(1)>0)
        t_sec(1) = 0;
        if(t_sec(2)>0)
            t_sec(2) = 0;
        end
    end
    count = size(escapeIndicesMS, 1);
    t_sec = resample(t_sec, count, length(t_sec));
    epm_n_sec_(jj) = t_sec(end, 1);
    % Activity
    try
        zc_raw = zscore(neural_data.C_raw, [], 2);
    catch
        zc_raw = zscore(neural_data, [], 2);
    end
    % Total time in center and average
    % t_m = Tracking.PercentTimeCenterEPM*n_sec;
    epm_t_m_p_(jj) = Tracking.PercentTimeCenterEPM;
    
    % Total time in closed and % of time spent in closed
    % t_c = Tracking.PercentTimeClosedEPM*n_sec;
    epm_t_c_p_(jj) = Tracking.PercentTimeClosedEPM;
    
    % Total time in closed and % of time spent in closed
    % t_o = Tracking.PercentTimeOpenEPM*n_sec;
    epm_t_o_p_(jj) = Tracking.PercentTimeOpenEPM;
    
    epm_ca_count_(jj) = size(closedArmFrameMS, 1);
    epm_oa_count_(jj) = size(openArmFrameMS, 1);
    
    epm_ca_count_p_(jj) = epm_ca_count_(jj)/(epm_ca_count_(jj)+epm_oa_count_(jj));
    epm_oa_count_p_(jj) = epm_oa_count_(jj)/(epm_ca_count_(jj)+epm_oa_count_(jj));
    
    % Total activity for cell ensembles
    % neither
    zc_neie = zc_raw((celltypes(1:size(zc_raw, 1))=='n'&keepind(1:size(zc_raw, 1))==1), :);
    % open
    zc_opne = zc_raw((celltypes(1:size(zc_raw, 1))=='o'&keepind(1:size(zc_raw, 1))==1), :);
    % closed
    zc_csde = zc_raw((celltypes(1:size(zc_raw, 1))=='c'&keepind(1:size(zc_raw, 1))==1), :);
    
    epm_zc_neie_mean_(jj) =  mean(mean(zc_neie, 2));
    epm_zc_opne_mean_(jj) =  mean(mean(zc_opne, 2));
    epm_zc_csde_mean_(jj) =  mean(mean(zc_csde, 2));
    
    
    % Activity of cell ensembles in each arm
    epm_oa_zc_neie_mean_(jj) = mean(mean(zc_neie(:, openArmIndicesMS(1:size(zc_raw, 2))==1), 2));
    epm_oa_zc_opne_mean_(jj) = mean(mean(zc_opne(:, openArmIndicesMS(1:size(zc_raw, 2))==1), 2));
    epm_oa_zc_csde_mean_(jj) = mean(mean(zc_csde(:, openArmIndicesMS(1:size(zc_raw, 2))==1), 2));
    
    epm_ca_zc_neie_mean_(jj) = mean(mean(zc_neie(:, closedArmIndicesMS(1:size(zc_raw, 2))==1), 2));
    epm_ca_zc_opne_mean_(jj) = mean(mean(zc_opne(:, closedArmIndicesMS(1:size(zc_raw, 2))==1), 2));
    epm_ca_zc_csde_mean_(jj) = mean(mean(zc_csde(:, closedArmIndicesMS(1:size(zc_raw, 2))==1), 2));
    
    epm_frz_percent_(jj) = (sum(freezeIndicesMS))/size(freezeIndicesMS, 1);
    epm_str_percent_(jj) = (sum(stretchIndicesMS))/size(stretchIndicesMS, 1);
    if(epm_frz_percent_(jj)>0)
        epm_frz_count_(jj) = size(freezeFrameMS, 1);
    end
    if(epm_str_percent_(jj)>0)
        epm_str_count_(jj) = size(stretchFrameMS, 1);
    end
    %     if exist('Tracking.mouseVelMS', 'var')
    %         epm_vel_mean_(jj) = nanmean(Tracking.mouseVelMS);
    %     else
    dmp = (diff(Tracking.mouse_position, 1).^2);
    epm_vel_mean_(jj) = (nansum(sqrt(dmp(:, 1)+dmp(:, 2))))/epm_n_sec_(jj);
    %     end
    clearvars Tracking escapeIndicesMS stretchIndicesMS freezeIndicesMS
    clearvars zc_raw
    
    cd ..
    cd('Rat1');
    %if (jj == 6)
    %    cd('H10_M45_S10');
    %end
    
    load('neural_data.mat');
    load('Tracking.mat');
    load('BehaviorMS.mat');
    
    % Organize the data a bit:
    % Time in seconds
    timestamp = readtable('timestamp.dat');
    t_sec = table2array(timestamp(:, 'sysClock'));
    t_sec = (t_sec/1000);
    if(t_sec(size(timestamp, 1), 1) == 0)
        timestamp([size(timestamp, 1)], :) = [];
        t_sec = table2array(timestamp(:, 'sysClock'));
        t_sec = (t_sec/1000);
    end
    if(t_sec(1)>0)
        t_sec(1) = 0;
        if(t_sec(2)>0)
            t_sec(2) = 0;
        end
    end
    count = size(escapeIndicesMS, 1);
    t_sec = resample(t_sec, count, length(t_sec));
    rat_n_sec(jj) = t_sec(end, 1);
    n_sec = rat_n_sec(jj);
    % Frames
    frame = table2array(timestamp(:, 'frameNum'));
    frame = resample(frame, count, length(frame));
    % Activity
    try
        c_raw = neural_data.C_raw;
    catch
        c_raw = neural_data;
    end
    zc_raw = zscore(c_raw, [], 2);
    n_cel = size(c_raw, 1);
    
    % Get celltypes and keepind from epm to rat
    rat_celltypes = char(zeros(n_cel, 1));
    rat_keepind = zeros(n_cel, 1);
    for i = 1:length(cor_ref)
        ref_to_rat_i = find(cor_rat == cor_ref(i));
        if (~(isempty(ref_to_rat_i))&(ref_to_rat_i>0))
            rat_to_epm_i = cor_epm(ref_to_rat_i);
            if (rat_to_epm_i>0)&(rat_to_epm_i<size(celltypes, 1))&(ref_to_rat_i<n_cel)
                rat_celltypes(ref_to_rat_i) = celltypes(rat_to_epm_i);
                rat_keepind(ref_to_rat_i) = keepind(rat_to_epm_i);
            end
        end
    end
    rat_celltypes_all{jj} = rat_celltypes;
    
    % Get far and near
    dim_x = max(max(Tracking.mouse_positionMS(:, 1)), max(Tracking.rat_positionMS(:, 1)));
    dim_y = max(max(Tracking.mouse_positionMS(:, 2)), max(Tracking.rat_positionMS(:, 2)));
    dim_nea_r = sqrt((dim_x*dim_y*.2)/pi);
    dim_far = (dim_x*.2);
    % Make booleans for far and near
    nea = Tracking.distanceMouseRat<dim_nea_r;
    nea = logical(resample(double(nea), count, length(nea)));
    nea_all{jj} = nea;
    far = Tracking.mouse_positionMS(:, 1)<dim_far;
    far = logical(resample(double(far), count, length(far)));
    far_all{jj} = far;
    
    rat_nea_percent_(jj) = sum(nea)/count;
    rat_far_percent_(jj) = sum(far)/count;
    
    rat_nea_count_(jj) = sum(diff(nea) == 1);
    rat_far_count_(jj) = sum(diff(far) == 1);
    
    % Total activity for cell ensembles
    % neither
    zc_neie = zc_raw((rat_celltypes=='n'&rat_keepind==1), :);
    % open
    zc_opne = zc_raw((rat_celltypes=='o'&rat_keepind==1), :);
    % closed
    zc_csde = zc_raw((rat_celltypes=='c'&rat_keepind==1), :);
    
    rat_zc_neie_mean_(jj) =  mean(mean(zc_neie, 2));
    rat_zc_opne_mean_(jj) =  mean(mean(zc_opne, 2));
    rat_zc_csde_mean_(jj) =  mean(mean(zc_csde, 2));
    
    
    % Activity of cell ensembles in nea and far
    rat_nea_zc_neie_mean_(jj) = mean(mean(zc_neie(:, nea==1), 2));
    rat_nea_zc_opne_mean_(jj) = mean(mean(zc_opne(:, nea==1), 2));
    rat_nea_zc_csde_mean_(jj) = mean(mean(zc_csde(:, nea==1), 2));
    
    rat_far_zc_neie_mean_(jj) = mean(mean(zc_neie(:, far==1), 2));
    rat_far_zc_opne_mean_(jj) = mean(mean(zc_opne(:, far==1), 2));
    rat_far_zc_csde_mean_(jj) = mean(mean(zc_csde(:, far==1), 2));
    
    rat_frz_percent_(jj) = Tracking.freezeFrac;
    rat_apr_percent_(jj) = (sum(approachIndicesMS))/count;
    rat_str_percent_(jj) = (sum(stretchIndicesMS))/count;
    rat_esc_percent_(jj) = (sum(escapeIndicesMS))/count;
    if(rat_frz_percent_(jj)>0)
        rat_frz_count_(jj) = size(freezeFrameMS, 1);
    end
    if(rat_apr_percent_(jj)>0)
        rat_apr_count_(jj) = size(approachFrameMS, 1);
    end
    if(rat_str_percent_(jj)>0)
        rat_str_count_(jj) = size(stretchFrameMS, 1);
    end
    if(rat_esc_percent_(jj)>0)
        rat_esc_count_(jj) = size(escapeFrameMS, 1);
    end
end

rat_nea_count_percent_ = rat_nea_count_./(rat_nea_count_+rat_far_count_);
p_opne = n_opne./t_cells;
p_csde = n_csde./t_cells;
p_neie = 1-p_opne-p_csde;
ratio_opne_csde = n_opne./n_csde;

epm_ascr_csde_mean = [-0.3097, -0.2518, -0.3265, -0.2821, -0.3619, -0.2896, -0.2635];

epm_ascr_opne_mean = [0.5168, 0.2968, 0.3626, 0.3978, 0.3649, 0.3511, 0.3028];

epm_entries_o_to_c = [0.1429, 0.2300, 0.1404, 0.2230, 0.2151, 0.1370, 0.2000];

epm_entries_c_to_o = [0.1548, 0.2400, 0.1579, 0.2230, 0.2258, 0.1370, 0.2000];

epm_entries_o_to_o = [0.1071, 0.0300, 0.0965, 0.0288, 0.0215, 0.0274, 0.1143];

epm_entries_c_to_c = [0.5952, 0.5000, 0.6053, 0.5252, 0.5376, 0.6986, 0.4857];

all_info_table = table(transpose(sess), transpose(p_opne), transpose(p_csde), transpose(p_neie), transpose(ratio_opne_csde), ...
    ...
    transpose(epm_t_o_p_), transpose(epm_t_c_p_), ...
    transpose(epm_ascr_csde_mean), transpose(epm_ascr_opne_mean), ...
    transpose(epm_ca_count_), transpose(epm_oa_count_), transpose(epm_ca_count_p_), transpose(epm_oa_count_p_), ...
    transpose(epm_zc_opne_mean_), transpose(epm_zc_csde_mean_), transpose(epm_zc_neie_mean_),...
    transpose(epm_oa_zc_opne_mean_), transpose(epm_oa_zc_csde_mean_), 	transpose(epm_oa_zc_neie_mean_), ...
    transpose(epm_ca_zc_opne_mean_), transpose(epm_ca_zc_csde_mean_), 	transpose(epm_ca_zc_neie_mean_), ...
    transpose(epm_frz_percent_), transpose(epm_str_percent_), ...
    transpose(epm_frz_count_), transpose(epm_str_count_), ...
    transpose(epm_entries_c_to_o), transpose(epm_entries_o_to_c), transpose(epm_entries_c_to_c), transpose(epm_entries_o_to_o), transpose(epm_n_sec_), ...
    ...
    transpose(rat_zc_opne_mean_), transpose(rat_zc_csde_mean_), transpose(rat_zc_neie_mean_),...
    transpose(rat_nea_zc_opne_mean_), transpose(rat_nea_zc_csde_mean_), 	transpose(rat_nea_zc_neie_mean_), ...
    transpose(rat_far_zc_opne_mean_), transpose(rat_far_zc_csde_mean_), 	transpose(rat_far_zc_neie_mean_), ...
    rat_far_percent_, rat_far_count_, ...
    rat_nea_percent_, rat_nea_count_, rat_nea_count_percent_, ...
    transpose(rat_frz_percent_), transpose(rat_frz_count_), ...
    transpose(rat_esc_percent_), transpose(rat_esc_count_), ...
    transpose(rat_str_percent_), transpose(rat_str_count_), ...
    transpose(rat_apr_percent_), transpose(rat_apr_count_), transpose(rat_n_sec_), 'VariableNames', ...
    {'Animal', 'percentCellsOpne', 'percentCellsCsde', 'percentCellsNeie', 'ratioOpneOverCsde', ...
    ...
    'EPMtimeOAPercent', 'EPMtimeCAPercent' ...
    'EPMarmScoreCsde', 'EPMarmScoreOpne', ...
    'EPMCACount', 'EPMOACount', 'EPMCACountPercent', 'EPMOACountPercent', ...
    'EPMOpneActivity', 'EPMCsdeActivity', 'EPMNeieActivity', ...
    'EPMOpneOAActivity', 'EPMCsdeOAActivity', 'EPMNeieOAActivity', ...
    'EPMOpneCAActivity', 'EPMCsdeCAActivity', 'EPMNeieCAActivity', ...
    'EPMfrzPercent', 'EPMfrzCount', ...
    'EPMstrPercent', 'EPMstrCount', ...
    'EPMentriesCAtoOA', 'EPMentriesOAtoCA', 'EPMentriesCAtoCA', 'EPMentriesOAtoOA', 'EPMSessDuration'...
    ...
    'RATOpneActivity', 'RATCsdeActivity', 'RATNeieActivity', ...
    'RATOpneNearActivity', 'RATCsdeNearActivity', 'RATNeieNearActivity', ...
    'RATOpneFarActivity', 'RATCsdeFarActivity', 'RATNeieFarActivity', ...
    'RATfarPercent', 'RATfarCount', ...
    'RATneaPercent', 'RATneaCount', 'RATneaCountPercent', ...
    'RATfrzPercent', 'RATfrzCount', ...
    'RATescPercent', 'RATescCount', ...
    'RATstrPercent', 'RATstrCount', ...
    'RATaprPercent', 'RATaprCount', 'RATSessDuration'});