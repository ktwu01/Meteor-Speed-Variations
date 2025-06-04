





% load original data
clear all;clc;close all;
addpath '/glade/u/home/wukoutian/0Res/Functions'

% % Update meteor parameters
% directoryPath = 'D:\0lrn\0Res\Data\';
% baseFilename = 'A0_Stations_Para_';
% WKT_runLatestVersion(directoryPath, baseFilename);
load('newmetsp(20140315-20231217).mat');

% % % % config
% site = 'ALTMR' % 70.0◦N, 23.3◦E
% years = 2016:2020;

% site = 'MHMR'
% years = 2011:2021;

% site = 'BJMR'
% years = 2011:2023;

site = 'MCMR'
years = 2014:2023;

% site = 'WHMR'
% years = 2021:2022;

% site = 'KMMR'
% years = 2011:2014;
% % % % % % % % % % % Note: KM do not have range and rd!!!
%
% % Know station name, find station id
% % site:  station name, ST_id: station id in D:\0lrn\0Res\Data\MR_Stations_Para.mat
% ST_id = WKT_findStringMatch(site, STs, sites)
% % STs{ST_id}
% % sites{ST_id}
%
% DAY_PER = datenum(years(end),12,31)- datenum(years(1),1,1)+1;
%
% % speed contourf data of each day of N years
% loadfilename1 = ...
%     sprintf('D:\\0lrn\\0Res\\Data\\%s_sp_dis_daily_%d_%d.mat', STs{ST_id}, years(1), years(end))
% load(loadfilename1, 'mcspdisNyrs');
%
% loadfilename2 = ...
%     sprintf('D:\\0lrn\\0Res\\Data\\%s_MetData_%d_%d.mat', STs{ST_id}, years(1), years(end))
% load(loadfilename2, 'number');
%
% if(1)
%     % Indices where 'number' values are below 1k
%     low_number_indices = find(number < 2000);
%     low_number_indices = low_number_indices(:);
%     mcspdisNyrs(:, [low_number_indices]) = nan;
% end
%

newmetsp = newmetsp'; % not start from 2014-01-01
StartIDD0 = datenum(2014,3,15)- datenum(2014,1,1); % Start Offset
EndIDD0 = datenum(2023,12,31)- datenum(2023,12,17); % End Offset
mcspdisNyrs =  nan(70,3565+StartIDD0+EndIDD0);
% mcspdisNyrs(:,1:StartIDD0) = nan; % start from 2014-01-01
mcspdisNyrs(:,StartIDD0+1:end-EndIDD0) = newmetsp; % start from 2014-01-01
% mcspdisNyrs(:,StartIDD0+1:end-EndIDD0) = nan; % start from 2014-01-01

% % a little test
% figure;
% contourf(1:DAY_PER,0:2:80,mcspdisNyrs,100,'Linestyle','none');colorbar;
% ylim([10 70])


% % % %% Calendar data
% Calendar data. %% use ShCal to fill in this data, and then plot
% START, END, MAX Date
calendars = {
    2014, {'QUA', [362, 12, 3]; 'ETA', [109, 148, 126]; 'SDA', [193, 235, 211]; 'GEM', [338, 351, 348]};
    2015, {'QUA', [362, 12, 4]; 'ETA', [109, 148, 126]; 'SDA', [193, 235, 211]; 'GEM', [338, 351, 348]};
    2016, {'QUA', [362, 12, 4]; 'ETA', [109, 148, 125]; 'ARI', [134, 175, 158]; 'SDA', [193, 235, 211]; 'GEM', [339, 352, 349]};
    2017, {'QUA', [362, 12, 3]; 'ETA', [109, 148, 126]; 'ARI', [134, 175, 158]; 'SDA', [193, 235, 211]; 'GEM', [339, 352, 349]};
    2018, {'QUA', [362, 12, 3]; 'ETA', [109, 148, 126]; 'ARI', [134, 175, 158]; 'SDA', [193, 235, 211]; 'GEM', [338, 351, 348]};
    2019, {'QUA', [362, 12, 3]; 'ETA', [109, 148, 126]; 'ARI', [134, 175, 158]; 'SDA', [193, 235, 211]; 'GEM', [338, 351, 348]};
    2020, {'QUA', [362, 12, 4]; 'ETA', [109, 148, 125]; 'ARI', [134, 175, 158]; 'SDA', [193, 235, 210]; 'GEM', [339, 354, 349]};
    2021, {'QUA', [362, 12, 3]; 'ETA', [109, 148, 125]; 'ARI', [134, 175, 158]; 'SDA', [193, 235, 211]; 'GEM', [339, 354, 349]};
    2022, {'QUA', [362, 12, 3]; 'ETA', [109, 148, 126]; 'ARI', [134, 175, 158]; 'SDA', [193, 235, 211]; 'GEM', [339, 354, 349]};
    2023, {'QUA', [362, 12, 4]; 'ETA', [109, 148, 126]; 'ARI', [134, 175, 158]; 'SDA', [193, 235, 211]; 'GEM', [339, 354, 349]};
    };

% % % % % % % Added: DSX
% % % % % % calendars = {
% % % % % %     2014, {'QUA', [362, 12, 3]; 'ETA', [109, 148, 126]; 'SDA', [193, 235, 211]; 'GEM', [338, 351, 348]};
% % % % % %     2015, {'QUA', [362, 12, 4]; 'ETA', [109, 148, 126]; 'SDA', [193, 235, 211]; 'GEM', [338, 351, 348]};
% % % % % %     2016, {'QUA', [362, 12, 4]; 'ETA', [109, 148, 125]; 'ARI', [134, 175, 158]; 'SDA', [193, 235, 211]; 'DSX', [252, 282, 270]; 'GEM', [339, 352, 349]};
% % % % % %     2017, {'QUA', [362, 12, 3]; 'ETA', [109, 148, 126]; 'ARI', [134, 175, 158]; 'SDA', [193, 235, 211]; 'DSX', [252, 282, 270]; 'GEM', [339, 352, 349]};
% % % % % %     2018, {'QUA', [362, 12, 3]; 'ETA', [109, 148, 126]; 'ARI', [134, 175, 158]; 'SDA', [193, 235, 211]; 'DSX', [252, 282, 270]; 'GEM', [338, 351, 348]};
% % % % % %     2019, {'QUA', [362, 12, 3]; 'ETA', [109, 148, 126]; 'ARI', [134, 175, 158]; 'SDA', [193, 235, 211]; 'DSX', [252, 282, 271]; 'GEM', [338, 351, 348]};
% % % % % %     2020, {'QUA', [362, 12, 4]; 'ETA', [109, 148, 125]; 'ARI', [134, 175, 158]; 'SDA', [193, 235, 210]; 'DSX', [252, 282, 270]; 'GEM', [339, 354, 349]};
% % % % % %     2021, {'QUA', [362, 12, 3]; 'ETA', [109, 148, 125]; 'ARI', [134, 175, 158]; 'SDA', [193, 235, 211]; 'DSX', [252, 282, 270]; 'GEM', [339, 354, 349]};
% % % % % %     2022, {'QUA', [362, 12, 3]; 'ETA', [109, 148, 126]; 'ARI', [134, 175, 158]; 'SDA', [193, 235, 211]; 'DSX', [252, 282, 270]; 'GEM', [339, 354, 349]};
% % % % % %     2023, {'QUA', [362, 12, 4]; 'ETA', [109, 148, 126]; 'ARI', [134, 175, 158]; 'SDA', [193, 235, 211]; 'DSX', [252, 282, 270]; 'GEM', [339, 354, 349]};
% % % % % %     };

figure; hold all;
whitefig;
% figpos([1 0.5])
% figpos([1 1])
figpos([0.5 1])
% fs = 12; % fontsize
fs = 19; % fontsize, FastX

%-------------------------------------------------------
vert_gap = 0.01;       horz_gap = 0;
lower_marg = 0.05;     upper_marg = 0.03;
left_marg = 0.17;      right_marg = 0.12;

Plot_Start_Yr = years(1);
Plot_End_Yr = years(end);
Target_Years = Plot_Start_Yr: Plot_End_Yr;
rows = length(Target_Years); cols = 1;
subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

% for Target_Year = 2018
for Target_Year = Target_Years

    % if Target_Year = 2015

    % speed contourf of Nyrs: method 4
    Start_IDX = datenum(Target_Year,1,1)- datenum(years(1),1,1)+1;
    End_IDX = datenum(Target_Year,12,31)- datenum(years(1),1,1)+1;

    % smooth period
    % SmtPrd = 30
    % sm1yr = nan(41, End_IDX-Start_IDX+1);
    mcspAno1yr = mcspdisNyrs(:,Start_IDX:End_IDX);
    % number1yr = number(Start_IDX:End_IDX);
    %
    % SG = 1;
    %
    % if(SG ==0)
    %     for i = 1:41
    %         sm1yr(i,:) = smoothdata (mcspdis1yr (i,:),'movmedian',SmtPrd);
    %     end
    %
    % end
    % % % % Wrong: sm1yr = smoothdata(mcspdis1yr, 'movmedian', SmtPrd);
    % if(SG ==1) % bad method
    %     degree = 3;
    %     SGSmtPrd = SmtPrd+1
    %     for i = 1:41
    %         sm1yr(i,:) = sgolayfilt(mcspdis1yr (i,:), degree, SGSmtPrd);
    %     end
    % end
    %
    % mcspAno1yr = (mcspdis1yr-sm1yr)./sm1yr;
    %
    % if Target_Year == 2016||Target_Year == 2018
    %     if(1)
    %         % Indices where 'number' values are below 1k
    %         low_number_indices = find(number1yr < 1000);
    %         low_number_indices = low_number_indices(:);
    %         % mcspAno1yr(:, [low_number_indices-1, low_number_indices, low_number_indices+1]) = nan;
    %
    %         try
    %             mcspAno1yr (:,low_number_indices-1 ) = nan;
    %         catch
    %         end
    %         mcspAno1yr (:,low_number_indices ) = nan;
    %         try
    %             mcspAno1yr (:,low_number_indices+1 ) = nan;
    %         catch
    %         end
    %
    %         try
    %             if Target_Year~=2016 & Target_Year~=2020 % "run nian" year
    %                 mcspAno1yr(:, 366) = [];
    %             end
    %         catch
    %             % do nothing
    %         end
    %     end
    % end
    %
    % % clear Inf anomaly value
    % mcspAno1yr (mcspAno1yr==Inf) = nan;
    %
    % if Target_Year == 2021
    %     mcspAno1yr(:, 189:191) = nan;
    % end
    %
    % if Target_Year == 2023
    %     mcspAno1yr(:, 176:179) = nan;
    % end

    ax = Target_Year-Plot_Start_Yr+1;
    subplot(rows,cols,ax)

    clen=100; cmap = nph_saturate(cbrew('nph_RainbowWhite',clen),1.7);
    hold on; colormap(gca,cmap);

    % !!!!!!!!!!!!!!!!!!!!!!!!!!!
    WTF = 10:79

    contourf(0:End_IDX-Start_IDX,WTF,mcspAno1yr,100,'Linestyle','none');

    %%%% LIMITS and ticks
    ylim([20 70])
    ytick(0:20:80);
    yminortick(0:5:80);
    xlim([0 366]);

    % clim([0.5 1.0]);
    % clim([0 1.0]);
    % clim([0 1.0]);

    hold on; axx=gca;
    axx.XTick = datenum(0,1:12,1)-datenum(0,01,01);
    axx.XMinorTick = 'off';
    axx.XTickLabel = {};

    % datetick('x','m','keepticks','keeplimits')
    hold on; nph_text([-0.13-0.09 0.85],['(' alphabet(ax) ')'],'fontsize',1.5*fs);
    % hold on; nph_text([-0.12 0.50],['(' num2str(Target_Year) ')'],'fontsize',1.5*fs);
    hold on; nph_text([-0.13-0.09 0.50],num2str(Target_Year),'fontsize',1.5*fs);

    box on;

    %% calendar dark
    % ... [Your existing code before this block] ...
    calendar_year = cell2mat(calendars(:,1));
    idx = find(calendar_year == Target_Year);

    if isempty(idx)
        warning(['No meteor shower data found for year ', num2str(Target_Year)]);
        continue; % Skip to the next iteration of the loop
    end

    current_year_meteor_showers = calendars{idx, 2};

    % Plot meteor showers for the current year
    for m = 1:size(current_year_meteor_showers, 1)
        shower_name = current_year_meteor_showers{m, 1};
        date_range = current_year_meteor_showers{m, 2}(1:2);
        max_date = current_year_meteor_showers{m, 2}(3);

        % Convert Dec dates of the previous year to the current year's
        % start
        if date_range(1) > 360
            date_range(1) = date_range(1) - 365;
        end
        if date_range(2) > 360
            date_range(2) = date_range(2) - 365;
        end
        if max_date > 360
            max_date = max_date - 365;
        end

        % Create patch
        p = patch('XData', [date_range(1), date_range(2), date_range(2), date_range(1)], ...
            'YData', [min(gca().YLim), min(gca().YLim), max(gca().YLim), max(gca().YLim)], ...  % Assuming y-axis range of 0 to 80
            'FaceColor', 'k', ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.1);
        % 'FaceColor', [0, 0, 0], ...  % RGB for black
        % % % Add annotation for the maximum date of the meteor shower %
        % y_position = 75; % A value slightly above the maximum y-axis
        % value % text(max_date, y_position, shower_name,
        % 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        % % % Add quiver to indicate the peak of the meteor shower
        % % quiver(max_date, 10, 0, 20, 0, 'k', 'LineWidth', 2, ...
        % 'MaxHeadSize', 0.5);
        % 'k' denotes black color, adjusted line width and arrowhead size
        xline(max_date, 'k', 'LineWidth', 0.5);
    end



    if Target_Year == Target_Years(1)
        % title('Seasonal variation of daily meteor speed distribution anomaly')
        Select_idx  = 4;
        current_year_meteor_showers = calendars{Select_idx, 2};

        for m = 1:size(current_year_meteor_showers, 1)
            shower_name = current_year_meteor_showers{m, 1};
            date_range = current_year_meteor_showers{m, 2};

            % Handle December dates from the previous year
            date_range(date_range > 360) = date_range(date_range > 360) - 365;

            % Calculate the x and y positions for the text
            x_position = mean(date_range);
            y_position = 73; % A value slightly above the maximum y-axis value

            % Add the text
            text(x_position, y_position, shower_name, 'fontweight','bold', 'fontsize',fs, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end

    end

    if Target_Year == 2019
        %%%% LABELS
        ylabel('Speed (km/s)');

        %%%% COLORBAR
        hold on;
        cbar = nph_colorbar;
        cbar.FontSize = fs;
        cbar.TickDirection = 'out';
        Cunit = cbar.Position(3);
        cbar.Position = cbar.Position .* [1 1 0.8 1];
        cbar.Position = cbar.Position + [0.3*Cunit -5*Cunit 0 2*5*Cunit];
        cbar.Label.Rotation = 90;
        cbar.Label.String = ['Meteor number'];
        % cbar.Label.String = ['Relative variability'];
        % cbar.Label.String = ['Count rate'];
        % cbar.Label.String = ['\bf{' clabs{ax} '}'];
        cbar.Ticks = [0:50:300];
    end

    if Target_Year == 2023
        % months as ticks using text at 15th day
        ytix1 = min(gca().YLim)-(max(gca().YLim)-min(gca().YLim))*0.01;
        xtixMt = datenum(0,1:12,15)-datenum(0,01,01);
        for xt = xtixMt
            if inrange(xt,xlim)
                mn = monthname(month(xt),'mmm');
                hold on; text(xt,ytix1,mn,'fontsize',fs,'horizontalalignment','center','VerticalAlignment','top');
            end
        end
    end
end

for ax = 1:rows*cols
    subplot(rows,cols,ax);
    axx = gca;
    setfont(fs);
    set(gca,'linewi',0.5,'tickdir','out');
end

error

%% EXPORT FIG ==============================================================
set(gcf, 'Renderer', 'painters');
% figure_name = ['Figure11_v5'+string(Target_Year)];
figure_name = ['Fig11'];
print(gcf,figure_name,'-djpeg','-r600');
% saveas(gcf,figure_name,'svg');
disp('Figure Saved.')