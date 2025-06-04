




%%% WKTsp_con_yrs_2p_V0

%% YOU MAY CHANGE THE FIG RATIO BY ANOTHER SCREEN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ AND TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % num: start from 2014-01-01
% % % % % % est: start from 2014-01-01 0:00 UT

% % % % % % min(mcest)=1.7575e+03 (73 days * 24 hours=1752)

% % % % % % recorded meteors: start from 74rd day 5.5 hour of 2014

% load original data
clear all;clc;close all; time_sys1=clock;
% load('../../Data/McMRdata_2014_2023.mat','doy','speed');

load('../../Data/Mc_sp_dis_daily_2014_2023.mat');
load('../../Data/Mc_sp_dis_ave_2014_2023.mat');

load('../../Data/Mc_sp_GS2_daily_2014_2023.mat');
load('../../Data/Mc_sp_GS2_ave_2014_2023.mat');

addpath '../../Functions'

DAY_PER = datenum(2023,12,31)-datenum(2014,1,1)+1;

% %% calculate the monthly peak and width

% mtdts = [31 28 31 30 31 30 31 31 30 31 30 31];
% indx = 1;
% for mt = 1:12
%     % DO NOT use: find (doy = dt)
%     mtcmpinds = find(doy>indx-1 & doy<indx + mtdts(mt)- 1+0.5);
%     indx = indx + mtdts(mt) - 1;
% %     figure;
%     H=histogram(speed(mtcmpinds)./1000,0:0.1:80,'visible','off');
%     title(['The',num2str(mt),'month'])
%     %%%%%%%%%%%%%%%%%%%%%%%%%% fit
%     % step1: normalize the number of meteors observed in range of mtcmpinds
%     % so that we can control the parameters in an acceptable region
%     y = H.BinCounts/max(H.BinCounts);
%     x = H.BinEdges(1:end-1);
% %     bar(x,y);
% %     xlabel('Speed (km/s)')
% %     ylabel('N/N_{max}')
% %     
% %     ylim([0 1.1])
% %     ytick(0:0.1:1)
% %     yminortick(0:0.05:1)
% 
%     % A1, b1, c1, A2, b2, c2
%     f = fit(x.',y.','gauss2', ...% a1, mu1, sig1, a2, mu2, sig2
%                 'Lower',[1, 26.0, 10, 0.20, 51, 13], ...
%                 'Start',[1, 28.0, 12, 0.33, 54, 15], ... 
%                 'Upper',[1, 30.0, 14, 0.46, 57, 17]);
% 
% %     yi= f.a1*exp(-((x-f.b1)/f.c1).^2) + f.a2*exp(-((x-f.b2)/f.c2).^2);
% %     hold on;plot(x,yi,'linest','--','color','b','linewi',2);
% 
%     spa2mt(mt)=f.a2;
%     spa1mt(mt)=f.a1; 
%     spwt2mt(mt)=f.c2; 
%     sppk2mt(mt)=f.b2;
%     spwt1mt(mt)=f.c1;
%     sppk1mt(mt)=f.b1;
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTY PLOT!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; 
whitefig;
figpos([1 0.4]) 

%-------------------------------------------------------
vert_gap = 0.055;       horz_gap = 0.06;
lower_marg = 0.15;      upper_marg = 0.07;
left_marg = 0.047;      right_marg = 0.056;

rows = 1; cols = 45;
subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 22; % fontsize
sz = 10; % size of the dots

axs = {...
    [1:cols-6],...
    [cols-6:cols]...
    };

% %% barcolor bar color
% % Example: barcolor = '#0065a1';
% barcolor = '#009BCA';
% dashcolor= '#dfdfdf';
% sphistcolor1	 = '	#4B6099	';
% sphistcolor2	 = '	#5E71A9	';
% sphistcolor3	 = '	#8086BA	';
% sphistcolor4	 = '	#9692C5	';
% sphistcolor5	 = '	#B6A7D2	';
% sphistcolor6	 = '	#C7B3CC	';
% sphistcolor7	 = '	#D6B7C7	';
% sphistcolor8	 = '	#E4BEBB	';
% sphistcolor9	 = '	#F1C5A8	';
% sphistcolor10	 = '	#FAC897	';
% sphistcolor11	 = '	#EFB285	';
% sphistcolor12	 = '	#A8736B	';


%% colormap
clen=100;cmap = nph_saturate(cbrew('nph_alt_jet',clen),1.0);

ytickrate1 =  0.0399;
ytickrate2 =  0.1111;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Speed against years %% Speed Vs time (Meteors Number against Year)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(rows,cols,axs{1})
for i = 1:DAY_PER
    mcspdisNyrs(:,i)=mcspdisNyrs(:,i)/max(mcspdisNyrs(:,i));
end

mcspdisNyrs=mcspdisNyrs/max(max(mcspdisNyrs));

hold on; colormap(gca,cmap);
hold on; contourf(1:DAY_PER,0:2:80,mcspdisNyrs,100,'Linestyle','none');
clim([0 1])

hold on; S2 = scatter(1:DAY_PER,sppk2, sz, 'filled', 'MarkerFaceColor', 'r', 'LineWidth',0.2);% xlim([0 DAY_PER]);ylabel('Speed Peak_2');
hold on; S1 = scatter(1:DAY_PER,sppk1, sz, 'filled', 'MarkerFaceColor', 'k', 'LineWidth',0.2);% ,'DisplayName','First Peak');

% %% colorbar
% hold on; cbar=colorbar; %clim([0 1]);
% % cbar=nph_colorbar;
% % cbar.Label.String = ['Daily meteor number'];
% % cbar.Label.String = ['N/N_{max}'];
% cbar.Label.String = ['Daily normalized number'];
% cbar.FontSize = fs;
% % cbar.Position(3) = 0.4*cbar.Position(3); % width
% % cbar.Position(4) = 4.2*cbar.Position(4); %height
% % % cbar.Position(1) = cbar.Position(1)-cbar.Position(3); %x axs
% % cbar.Position(2) = cbar.Position(2)-2.6*cbar.Position(3);  %y axs

%%%% LIMITS and ticks
ylim([10 80])
% ylim([0 80])
ytick(0:10:80);
yminortick(0:2:80);

%%%% LABELS
years = {};
years = 2014:2023;

XLIM_DAY = datenum((max(years)),07,03)-datenum(min(years),01,01);
xlim([0 XLIM_DAY])

xtickYrs(gca, years, 'k');
xtickYrslabel(gca, years,fs);

ylabel('Speed (km/s)')
title ('Annual')

% xlim([0 3500]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Speed against years %% Speed Vs time (Meteors Number against Year)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(rows,cols,axs{2})

for i = 1:366
    mcspdisave(:,i)=mcspdisave(:,i)/max(mcspdisave(:,i));
end

%% contourf

hold on; colormap(gca,cmap);
hold on; contourf(0:365,0:2:80,mcspdisave,60,'Linestyle','none');

hold on; scatter(0:365,sppk1ave,sz,'filled','MarkerFaceColor','k');
hold on; scatter(0:365,sppk2ave,sz,'filled','MarkerFaceColor','r');

%%%% LIMITS and ticks
ylim([10 80])
ytick(0:10:80);
yminortick(0:2:80);
xlim([0 365]);

% grid on;
hold on; axx=gca;
axx.XTick = datenum(2014,1:12,1)-datenum(2014,01,01);
datetick('x','m','keepticks','keeplimits')

axx.XTickLabel = {};
% months as ticks using text at 15th day
ytix1 = min(gca().YLim)-(max(gca().YLim)-min(gca().YLim))*ytickrate1;
xtixMt = datenum(2014,[1 4 7 10 13],5)-datenum(2014,01,01);
for xt = xtixMt
    if inrange(xt,xlim)
        mn = monthname(month(xt),'mmm');
    hold on; text(xt,ytix1,mn,'fontsize',fs,'horizontalalignment','center','VerticalAlignment','top');
    end
end

ytix1 = min(gca().YLim)-(max(gca().YLim)-min(gca().YLim))*ytickrate1;
xtixMt = datenum(0,12,30);
for xt = xtixMt
    hold on; text(xt,ytix1,'Jan','fontsize',fs,'horizontalalignment','center','VerticalAlignment','top');
end

axx.XMinorTick = 'off';

% Year as lower numbers:
ytix2 = min(gca().YLim)-(max(gca().YLim)-min(gca().YLim))*ytickrate2;
xtix2 = datenum(2014,07,01)-datenum(2014,01,01);
yrstr = 'Years average';
text(xtix2,ytix2,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');

% ylabel('Speed (km/s)')
title('Seasonal')

% axx.YTick = '';
% axx.YMinorTick = 'on';

%%%% COLORBAR
hold on; 
clim([0 1]);
cbar = nph_colorbar;
cbar.Label.String = ['Normalized number'];
% cbar.Ticks = [-7 -5 -3 -1 0];
% cbar.Ticks = -8:2:0;
cbar.TickDirection = 'out';
cbar.Position = cbar.Position .* [1 1 0.8 1];
cbar.Position = cbar.Position + [0.5*cbar.Position(3) 0 0 0];
%             cbar.Label.String = ['\bf{' clabs{ax} '}'];
cbar.Label.Rotation = 90;
% box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GLOBAL FORMATTING FOR ALL AXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ax = 1:2
    subplot(rows,cols,axs{ax});
    axx = gca;

    setfont(fs);
    set(gca,'linewi',1,'tickdir','out');

    switch ax
        case {1}
            hold on; nph_text([-0.01 0.899],['(' alphabet(ax) ')'],'fontsize',1.5*fs); % ,'textborder','w');
        case 2
            hold on; nph_text([ 0.1 0.899],['(' alphabet(ax) ')'],'fontsize',1.5*fs); % ,'textborder','w');
    end
end
time_sys2=clock;
timerun_21=etime(time_sys2,time_sys1);
time21=num2str(floor(timerun_21));


error

%% EXPORT FIG ==============================================================
% Set renderer to "painters"
set(gcf, 'Renderer', 'painters');
% print(gcf,['WKTWKTsp_con_yrs_2p_V0_t=',time21,'s'],'-dpng','-r600');
% saveas(gcf,'Figure9','epsc')
% saveas(gcf,['WKTWKTsp_con_yrs_2p_V0_t=',time21,'s'],'svg')
figure_name = ['Fig9'];
print(gcf,figure_name,'-djpeg','-r600');

disp('Figure Saved.')


function xtickYrs(gca, years, dashcolor)

hold on;grid on; axx=gca;

% Just show the Januarys as major ticks
axx.XTick = datenum(min(years):(max(years)+1),01,01)-datenum(min(years),01,01);
axx.XMinorTick = 'on';
axx.XAxis.MinorTickValues = datenum(min(years),1:((length(years)+1)*12),01)-datenum(min(years),01,01);
datetick('x','m','keepticks','keeplimits');
axx.XTickLabel = {};

% add dashed lines to the Januarys except the first Jan
dashday = datenum((min(years)):(max(years)+1),01,01)-datenum(min(years),01,01);
for xt = dashday
    if inrange(xt,xlim)% && ~any(xt == axx.XTick)
    hold on; plot([xt xt],axx.YLim,'linest','--','color',dashcolor,'linewi',axx.LineWidth);
%     text(xt,-500./1000,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
    end
end

end



function xtickYrslabel(gca, years, fs)

ytickrate1 =  0.0399;
ytickrate2 =  0.1111;

hold on;grid on; axx=gca;

% January & Apr & Jul as labels using text:
ytix1 = min(gca().YLim)-(max(gca().YLim)-min(gca().YLim))*ytickrate1;
xtixJan = datenum(years(1),1:12:(length(years)*12),2)-datenum(min(years)-1,12,31);
for xt = xtixJan
%     if inrange(xt,xlim) % && ~any(xt == axx.XTick)
        mn = monthname(month(xt),'mmm');
        text(xt,ytix1,mn,'fontsize',fs,'horizontalalignment','center','VerticalAlignment','top');
%     end
end

xtixJul = datenum(years(1),7:12:(length(years)*12+7),2)-datenum(min(years)-1,12,31);
for xt = xtixJul
    if inrange(xt,xlim) % && ~any(xt == axx.XTick)
        mn = monthname(month(xt),'mmm');
        text(xt,ytix1,mn,'fontsize',fs,'horizontalalignment','center','VerticalAlignment','top');
    end
end

% % other months as labels using text:
% xtix = datenum(years(1),1:(length(years)*12),15)-datenum(min(years)-1,01,01);
% for xt = xtix
%     if inrange(xt,xlim) && ~any(xt == axx.XTick)
%         mn = monthname(month(xt));
%         text(xt,-500./1000,mn(1),'fontsize',fs,'horizontalalignment','center','VerticalAlignment','top');
%     end
% end

% Years as lower numbers:
ytix2 = min(gca().YLim)-(max(gca().YLim)-min(gca().YLim))*ytickrate2;
xtix2 = datenum(years,07,01)-datenum(min(years),01,01);
for xt = xtix2
    if inrange(xt,xlim) && ~any(xt == axx.XTick)
        yrstr = datestr(xt,'yyyy');
        yrstr = num2str(str2double(yrstr) + min(years));
%         bold number
%         text(xt,ytix2,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
        text(xt,ytix2,yrstr,'fontsize',1*fs,'horizontalalignment','center','VerticalAlignment','top');
    end
end

end

