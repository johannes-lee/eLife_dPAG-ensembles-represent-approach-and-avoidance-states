%% Example arbitrary correlation plot
% SMP, updated 07/23/20

% Time spent near rat (x axis) vs mean CLOSED activity NEAR
%     subplot(2, 2, 4);
col_cck = [85/255 45/255 134/255];
col_syn = [0/255 138/255 179/255];
figure('Name', 'OF signal and velocity correlations', 'Color', 'w', 'NumberTitle', 'off');
subplot(1, 2, 1);
t_x2 = transpose(sig_mean(:, 1:3));
a_y2 = transpose(vel_mean(:, 1:3));
t_x = t_x2;
a_y = a_y2;
c = col_cck; ...linspace(0, 10, size(t_x, 1));

t_ac_corrc = scatter(t_x, a_y, [], c, 'd', 'filled');
%   hold on
hold on
coeffs = polyfit(t_x, (a_y), 1);
% Get fitted values
fitX = [min(t_x) (max(t_x))];
fitY = polyval(coeffs, fitX);
hold on
line(fitX, fitY, 'LineStyle', '--', 'Color', [0.25, 0.25, 0.25]);
hold off
hold off
[r2, p2] = corr(t_x, (a_y),'Type','Spearman');
[r, p] = corr(t_x, (a_y),'Type','Pearson');
title({'\fontsize{11} {\bfSignal and Speed CCK}', ['\fontsize{9} r_{spea} = ' num2str(r2) ', p_{spea} = ' num2str(p2), ...
    ', ', 'r_{pear} = ' num2str(r) ', p_{pear} = ' num2str(p) ]...
    }, 'Fontweight', 'normal')
ylabel('mean speed');
xlabel('mean dF/F');
% xlim([0 ceil(max(t_x))]);
% ylim([(min(t_x)) (max(a_y))]);

subplot(1, 2, 2);
t_x2 = transpose(sig_mean(:, 5:end));
a_y2 = transpose(vel_mean(:, 5:end));
t_x = t_x2;
a_y = a_y2;
c = col_syn; ...linspace(0, 10, size(t_x, 1));
t_ac_corrc = scatter(t_x, a_y, [], c, 'd', 'filled');
%   hold on
hold on
coeffs = polyfit(t_x, (a_y), 1);
% Get fitted values
fitX = [min(t_x) (max(t_x))];
fitY = polyval(coeffs, fitX);
hold on
line(fitX, fitY, 'LineStyle', '--', 'Color', [0.25, 0.25, 0.25]);
hold off
hold off
[r2, p2] = corr(t_x, (a_y),'Type','Spearman');
[r, p] = corr(t_x, (a_y),'Type','Pearson');
title({'\fontsize{11} {\bfSignal and Speed Syn}', ['\fontsize{9} r_{spea} = ' num2str(r2) ', p_{spea} = ' num2str(p2), ...
    ', ', 'r_{pear} = ' num2str(r) ', p_{pear} = ' num2str(p) ]...
    }, 'Fontweight', 'normal')
ylabel('mean speed');
xlabel('mean dF/F');
% xlim([0 ceil(max(t_x))]);
% ylim([(min(t_x)) (max(a_y))]);

% subplot(2, 2, 3);
% t_x2 = transpose(time_ca);
% a_y2 = transpose(zsc_sig_oa);
% t_x = t_x2;
% a_y = a_y2;
% c = linspace(0, 10, size(t_x, 1));
% t_ac_corrc = scatter(t_x, a_y, [], c, 'd', 'filled');
% %   hold on
% hold on
% coeffs = polyfit(t_x, (a_y), 1);
% % Get fitted values
% fitX = [0 ceil(max(t_x))];
% fitY = polyval(coeffs, fitX);
% hold on
% line(fitX, fitY, 'LineStyle', '--', 'Color', [0.25, 0.25, 0.25]);
% hold off
% hold off
% [r2, p2] = corr(t_x, (a_y),'Type','Spearman');
% [r, p] = corr(t_x, (a_y),'Type','Pearson');
% title({'\fontsize{11} {\bfOA activity, CA exploration}', ['\fontsize{9} r_{spea} = ' num2str(r2) ', p_{spea} = ' num2str(p2), ...
%     ', ', 'r_{pear} = ' num2str(r) ', p_{pear} = ' num2str(p) ]}, 'Fontweight', 'normal')
% ylabel('z-scored dF/F in OA');
% xlabel('% time in CA');
% xlim([0 ceil(max(t_x))]);
% % ylim([(min(t_x)) (max(a_y))]);
% 
% subplot(2, 2, 4);
% t_x2 = transpose(time_oa);
% a_y2 = transpose(zsc_sig_ca);
% t_x = t_x2;
% a_y = a_y2;
% c = linspace(0, 10, size(t_x, 1));
% t_ac_corrc = scatter(t_x, a_y, [], c, 'd', 'filled');
% %   hold on
% hold on
% coeffs = polyfit(t_x, (a_y), 1);
% % Get fitted values
% fitX = [0 ceil(max(t_x))];
% fitY = polyval(coeffs, fitX);
% hold on
% line(fitX, fitY, 'LineStyle', '--', 'Color', [0.25, 0.25, 0.25]);
% hold off
% hold off
% [r2, p2] = corr(t_x, (a_y),'Type','Spearman');
% [r, p] = corr(t_x, (a_y),'Type','Pearson');
% title({'\fontsize{11} {\bfCA activity, OA exploration}', ['\fontsize{9} r_{spea} = ' num2str(r2) ', p_{spea} = ' num2str(p2), ...
%     ', ', 'r_{pear} = ' num2str(r) ', p_{pear} = ' num2str(p) ]}, 'Fontweight', 'normal')
% ylabel('z-scored dF/F in CA');
% xlabel('% time in OA');
% xlim([0 ceil(max(t_x))]);
