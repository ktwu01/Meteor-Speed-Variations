% %%





% %% WKT_sp_diff_mtd4_V1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ AND TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % num: start from 2014-01-01
% % % % % % est: start from 2014-01-01 0:00 UT

% % % % % % 24*(datenum(2022,12,31)- datenum(2014,01,01)+1)=max(mcest) 
% % % % % % min(mcest)=1.7575e+03 (73 days * 24 hours=1752)

% % % % % % recorded meteors: start from 74rd day 5.5 hour of 2014

% load original data
clear all;clc;close all; time_sys1=clock;

load('../../Data/Mc_sp_dis_ave_2014_2023.mat');

addpath '../../Functions'

% DAY_PER = datenum(2023,12,31)-datenum(2014,1,1)+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTY PLOT!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; 
whitefig;

figpos([1 0.5])

%-------------------------------------------------------
vert_gap = 0.0;        horz_gap = 0;
lower_marg = 0.06;     upper_marg = 0.03;
left_marg = 0.05;      right_marg = 0.08;

rows = 1; cols = 1;
subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------
fs = 30; % fontsize
% sz = 20; % size of the dots

% % axs = {1:cols cols+1:cols*2 cols*2+1:cols*3}; % the same widths
% axs = {[1:cols cols*1+1:cols*2] [cols*2+1:cols*3 cols*4+1:cols*5]};

%% calculate

%% smooth period
SmtPrd = 30

for i = 1:41
    smave(i,:) = smoothdata (mcspdisave(i,:),'movmedian',SmtPrd);
end

mcspdisave = (mcspdisave-smave)./smave;
%     subplot 211;
%     hold on;plot(1:366,msc0);
%     hold on;plot(1:366,mscfit);
%     title([num2str(i*2-2) 'km/s'])
% 
%     subplot 212;hold on;plot(1:366,remsc(i,:));
%     title([num2str(i*2-2) 'km/s'])

%% smooth by rows
if(0)
% figure;
for i=1:366
%     mscc0=remsc(:,i);
    mscc0=mcspdisave(:,i);
    SmtPrd = 30; % smooth period
    mscfit=smooth(mscc0,SmtPrd);
    remscc(:,i)=mscc0-0.5*mscfit;
end
    mcspdisave =  remscc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Speed anomaly against 9 years %% Speed Vs time (Meteors Number against Year)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subplot(rows,cols,axs{1})
subplot(rows,cols,1)

clen=100; cmap = nph_saturate(cbrew('nph_RainbowWhite',clen),1);
% clen=100; cmap = nph_saturate(cbrew('nph_BuOrRd',clen),0.7); % can not be CyclicEdge!!
% cmap = jet;
hold on; colormap(gca,cmap);
% hold on; contourf(0:365,0:2:80,(remsc),100,'Linestyle','none');
% hold on; contourf(0:365,20:2:70,mcspdisave(10:1:35,:),100,'Linestyle','none');
hold on; contourf(0:365,0:2:80,mcspdisave,100,'Linestyle','none');

% %%% LIMITS and ticks
ylim([20 70])
ytick(0:10:80);
yminortick(0:2:80);
xlim([0 366]);

% %%% LABELS
ylabel('Speed (km/s)');
% title('Seasonal variation of daily meteor speed distribution anomaly')

% %%% COLORBAR
hold on;
cbar = nph_colorbar;
% cbar.Label.String = ['Normalized number'];
cbar.Label.String = ['Relative variability'];

clim([0 1.0]);
% clim([0 0.5]);
cbar.Ticks = [0:0.1:1.0];
% cbar.Ticks = -8:2:0;
cbar.TickDirection = 'out';
cbar.Position = cbar.Position .* [1 1 0.8 1];
cbar.Position = cbar.Position + [0.5*cbar.Position(3) 0 0 0];
%             cbar.Label.String = ['\bf{' clabs{ax} '}'];
cbar.Label.Rotation = 90;
box on;

hold on; axx=gca;
axx.XTick = datenum(2014,1:12,1)-datenum(2014,01,01);
datetick('x','m','keepticks','keeplimits')

axx.XTickLabel = {};

% months as ticks using text at 15th day
ytix1 = min(gca().YLim)-(max(gca().YLim)-min(gca().YLim))*0.01;
xtixMt = datenum(2014,1:12,15)-datenum(2014,01,01);
for xt = xtixMt
    if inrange(xt,xlim)
        mn = monthname(month(xt),'mmm');
    hold on; text(xt,ytix1,mn,'fontsize',fs,'horizontalalignment','center','VerticalAlignment','top');
    end
end

axx.XMinorTick = 'off';

%% legend # 1: numbers with circles indicate meteor showers
xshw = 2; yshw = 41+4; % nshw = {'1'};

swsz = fs*65;
scatter(xshw+2,yshw+9,swsz,'MarkerEdgeColor','k','MarkerFaceColor','none','LineWidth',2);
text(xshw,yshw+9,'1','FontSize',fs*1.4,'FontWeight','bold');
% text(xshw-4,yshw+5,"|",'FontSize',fs*1.4,'FontWeight','bold');
text(xshw-2,yshw+6,"\downarrow",'FontSize',fs*1.4,'FontWeight','bold');

%% legend # 2: define meteor showers

xshw(1) =  datenum('6-May-2014','dd-mmm-yyyy')- datenum(2014,01,01);
yshw(1) = 66-18; % ETA
% nshw{1} = {'2'};

scatter(xshw,yshw,swsz,'MarkerEdgeColor','k','MarkerFaceColor','none','LineWidth',2);
text(xshw-2,yshw,'2','FontSize',fs*1.4,'FontWeight','bold');
% text(xshw-4,yshw+5,"|",'FontSize',fs*1.4,'FontWeight','bold');
text(xshw-3,yshw+4,"\uparrow",'FontSize',fs*1.4,'FontWeight','bold');

% xshw(3) =  datenum('21-Oct-2014','dd-mmm-yyyy')- datenum(2014,01,01)-1;
% yshw(3) = 66; % ORI
% nshw{3} = {'Z'};

%% legend # 3-5 : numbers with circles indicate meteor showers

xshw(1) =  datenum('7-Jun-2014','dd-mmm-yyyy')- datenum(2014,01,01)+4;
yshw(1) = 38+4;
nshw{1} = {'3'};
% nshw{3} = {'ARI',' Max: June 8'};

xshw(2) =  datenum('30-Jul-2014','dd-mmm-yyyy')- datenum(2014,01,01)-3;
yshw(2) = 41+4;
nshw{2} = {'4'};
% nshw{5} = {'SDA',' Max: Jul 30'};

xshw(3) = datenum('14-Dec-2014','dd-mmm-yyyy')- datenum(2014,01,01)-1;
yshw(3) = 35+4;
nshw{3} = {'5'};
% nshw{9} = {'GEM',' Max: Dec 14'};


%% legend: numbers with circles indicate meteor showers

% text(xshw-4,yshw+5,"|",'FontSize',fs*1.4,'FontWeight','bold');
text(xshw-2,yshw+6,"\downarrow",'FontSize',fs*1.4,'FontWeight','bold');

scatter(xshw+1,yshw+9,swsz,'MarkerEdgeColor','k','MarkerFaceColor','none','LineWidth',2);
text(xshw-1,yshw+9,nshw,'FontSize',fs*1.4,'FontWeight','bold');

% for xt = xshw
%     if inrange(xt,xlim)% && ~any(xt == axx.XTick)
% %% legend: lines indicate meteor showers
%     hold on; plot([xt xt],axx.YLim,'linest','--','color','k','linewi',axx.LineWidth);
%     end
% end

%% scatter speed peaks
if (0)
hold on; scatter(0:365,sppk1ave,sz,'filled','MarkerFaceColor','k');
hold on; scatter(0:365,sppk2ave,sz,'filled','MarkerFaceColor','r');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GLOBAL FORMATTING FOR ALL AXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ax = 1
%     subplot(rows,cols,axs{ax});
    subplot(rows,cols,ax);
    axx = gca;
    setfont(fs);
    set(gca,'linewi',1,'tickdir','out');
%     hold on; nph_text([-0.125 0.855],['(' alphabet(ax) ')'],'fontsize',1.5*fs);
%     hold on; nph_text([0.025 0.865],['(' alphabet(ax) ')'],'fontsize',1.5*fs,'textborder','w');
end

time_sys2=clock;
timerun_21=etime(time_sys2,time_sys1);
time21=num2str(floor(timerun_21));

error
%% EXPORT FIG ==============================================================
% Set renderer to "painters"
set(gcf, 'Renderer', 'painters');
figure_name = ['Fig10'];
print(gcf,figure_name,'-djpeg','-r600');
% print(gcf,['WKT_sp_diff_mtd4_V0_t=',time21,'s'],'-dpng','-r600')
% saveas(gcf,'Figure10','epsc')
% saveas(gcf,['WKT_sp_diff_mtd4_V0_t=',time21,'s'],'svg')
disp('Figure Saved.')
