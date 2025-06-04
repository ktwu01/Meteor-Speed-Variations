



%%% azi_lt_surf_con_1_4p_V1_4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ AND TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % num: start from 2014-01-01
% % % % % % est: start from 2014-01-01 0:00 UT

% % % % % % 24*(datenum(2022,12,31)- datenum(2014,01,01)+1)=max(mcest) 
% % % % % % min(mcest)=1.7575e+03 (73 days * 24 hours=1752)

% % % % % % recorded meteors: start from 74rd day 5.5 hour of 2014

clear all; clc; close all; time_sys1=clock;
% load original data
load('D:\0lrn\00Res\Data\OptionalMMRdata2014_2022.mat', ...
    'azimuth','doy','localtime','mcest','speed');
load(['D:\0lrn\00Res\Data\ids2014_2022.mat']);

addpath 'D:\0lrn\00Res\Functions'

site = 'MC'; % Mengcheng

% day range
years = 2014:2022;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTY PLOT!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------
figure; hold all;
whitefig;
% figpos([0.75 0.8])+ax(5)
%-------------------------------------------------------
rows = 3; cols = 8;
figpos([1 1]);
vert_gap = 0.075*sqrt(5/rows); horz_gap = 0.05*sqrt(10/cols);
lower_marg = 0.05*sqrt(5/rows); upper_marg = 0.04*sqrt(5/rows);
left_marg = 0.045*sqrt(10/cols); right_marg = 0.035*sqrt(10/cols);
subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);
%--------------------------------------------------------
fs = 10; % fontsize
sz = 6; % scattersize
mks = 3; % markersize of meteor map

axs = {
    [1:8 9:16],
    17:18,
    19:20,
    21:22,
    23:24
    };

subplot(rows,cols,axs{1});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calcu contourf
%% inds: find 3 days centered on 2022,10,02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2022,10,02 localtime = 
selected_day = datenum(2022,10,02)-1;

daystrbefore = datestr((datenum(2022,10,02)-1),'mmm dd');
daystron = datestr((datenum(2022,10,02)),'mmm dd');
daystrafter = datestr((datenum(2022,10,02)+1),'mmm dd');

localtimelong = mcest+8; % local est time
inds = find(localtimelong>24*(selected_day- datenum(2014,01,01)) ...
    & localtimelong<24*(selected_day- datenum(2014,01,01)+3));
localtime1=localtimelong(inds)-24*(selected_day- datenum(2014,01,01));
%%%%%%%%%%%%  do not use: localtime1=localtime(inds); because 0<localtime<24
azimuth1 = azimuth(inds);
speed1=speed(inds)/1000;

azidtdis = zeros(72+1,41);
n=1;
for lt=0:1:72
    m=1;
    for az=0:9:360
        azidtdis(n,m)=length(find(azimuth1>az-4.5 & ...
            azimuth1<az+4.5 ...
            & localtime1>lt-0.5 & localtime1<lt+0.5));
        m=m+1;
    end
    n=n+1;
end

% %  normalize each row to unit
for lt = (0:1:72)+1
    azidtdis(lt,:) = azidtdis(lt,:)/max(azidtdis(lt,:));
end

% azidtdis = azidtdis./repmat(sqrt(sum(azidtdis.^2,2)),1,size(azidtdis,2));

azidtdis(:)  = smooth(azidtdis(:), 3);
azidtdis=azidtdis/max(max(azidtdis));

% fig

viewazi = -70;
viewele = 50;
% view([Azimuth,Elevation])
view([viewazi,viewele])

% colormap
clen=100;cmap = nph_saturate(cbrew('nph_alt_jet',clen),1.0);

hold on; colormap(gca,cmap);
title(datestr((datenum(2022,10,02)),'mmm dd yyyy'));
[xx,yy] = meshgrid (0:9:360,0:1:72);
surf(xx,yy,azidtdis,'FaceAlpha',1,'LineStyle','none');

hold on;
plot3(0:9:360,repmat(0,1,41),azidtdis(1,:),LineWidth=2,Color='r');
hold on;
plot3(repmat(0,1,72+1),0:1:72,azidtdis(:,1),LineWidth=2,Color='r');

xlim([0 360]);
xtick(0:60:360)
xminortick(0:10:360)

ylim([0 72])
ytick(0:12:72)
yticklabel({daystrbefore,'12',daystron,'12',daystrafter,'12',''})
yminortick(0:1:72)


zlim([0 1])
ztick(0:0.2:1)
zminortick(0:0.1:1)

%%%% LABELS
ylabel('Local time (hour)');
% zlabel('log_{10}(Hourly normalized number)');
zlabel('Hourly normalized number');
xlabel('Azimuth (degree)');

%%%% COLORBAR
% % % cbar = nph_colorbar;
% % % % cbar.Ticks = [-7 -5 -3 -1 0];
% % % % cbar.Ticks = -7.9:2:0;
% % % cbar.TickDirection = 'out';
% % % cbar.Position = cbar.Position .* [1 2 0.8 0.4];
% % % cbar.Position = cbar.Position + [-1.5*cbar.Position(3) -12*cbar.Position(3) 0 0];
% % % % % cbar.Label.String = ['\bf{' clabs{ax} '}'];
% % % cbar.Label.Rotation = 0;

grid on; % box on;

%     axx = gca;
%     setfont(fs);
%     set(gca,'linewi',0.5,'tickdir','out');
%     hold on; nph_text([-0.150 1.10],['(' alphabet(1) ')'],'fontsize',1.5*fs);

subplot(rows,cols,axs{2});azidtdis1=az_lt_contourf(id_pr,localtime, azimuth, 1);title('Jan 4');ylabel('Azimuth (degree)')
%     axx = gca;
%     setfont(fs);
%     set(gca,'linewi',0.5,'tickdir','out');
%     hold on; nph_text([-0.150 1.10],['(' alphabet(2) ')'],'fontsize',1.5*fs);

subplot(rows,cols,axs{3});azidtdis2=az_lt_contourf(id_se,localtime, azimuth, 2);title('Apr 4')
%     axx = gca;
%     setfont(fs);
%     set(gca,'linewi',0.5,'tickdir','out');
%     hold on; nph_text([-0.150 1.10],['(' alphabet(3) ')'],'fontsize',1.5*fs);
% 
subplot(rows,cols,axs{4});azidtdis3=az_lt_contourf(id_ap,localtime, azimuth, 3);title('Jul 4')
%     axx = gca;
%     setfont(fs);
%     set(gca,'linewi',0.5,'tickdir','out');
%     hold on; nph_text([-0.150 1.10],['(' alphabet(4) ')'],'fontsize',1.5*fs);

subplot(rows,cols,axs{5});azidtdis4=az_lt_contourf(id_ae,localtime, azimuth, 4);title('Oct 2')
%     axx = gca;
%     setfont(fs);
%     set(gca,'linewi',0.5,'tickdir','out');
%     hold on; nph_text([-0.150 1.10],['(' alphabet(5) ')'],'fontsize',1.5*fs);

cbar= nph_colorbar;
% cbar.Label.String = ['Hourly normalized number'];
    % nph_colorbar


for ax = 1:5
    subplot(rows,cols,axs{ax});
    axx = gca;
    setfont(fs);
    set(gca,'linewi',0.5,'tickdir','out');
    switch ax
        case 1
        hold on; nph_text([-0.015 0.87],['(' alphabet(ax) ')'],'fontsize',1.5*fs);
        case {2 3 4 5}
        hold on; nph_text([-0.17 1.1],['(' alphabet(ax) ')'],'fontsize',1.5*fs); %,'textborder','w');
    end
end



time_sys2=clock;
timerun_21=etime(time_sys2,time_sys1);



%% EXPORT FIG ==============================================================
% Set renderer to "painters"
set(gcf, 'Renderer', 'painters');
print(gcf,['WKTazi_lt_surf_con_1_4p_V1_4'],'-dpng','-r300');
% saveas(gcf,'WKTazi_lt_surf_con_1_4p_V1_4','epsc')
% saveas(gcf,'WKTazi_lt_surf_con_1_4p_V1_4','svg')
disp('Fig:WKTazi_lt_surf_con_1_4p_V1_4 Saved.')

function azidtdis=az_lt_contourf(id_ae,localtime, azimuth, px)
    % colormap
    clen=100;cmap = nph_saturate(cbrew('nph_alt_jet',clen),1.0);

    % calcu contourf
    localtime1=localtime(id_ae);
    azimuth1 = azimuth(id_ae);
    
    azidtdis = zeros(25,41);
    n=1;
    
    for lt=0:1:24
        m=1;
        for az=0:9:360
            azidtdis(n,m)=length(find(azimuth1>az-4.5 & ...
                azimuth1<az+4.5 ...
                & localtime1>lt-0.5 & localtime1<lt+0.5));
            m=m+1;
        end
        n=n+1;
    end
    
    % %  normalize each row to unit
    for lt = (0:1:24)+1
    azidtdis(lt,:) = azidtdis(lt,:)/max(azidtdis(lt,:));
    end

%     azidtdis = azidtdis./repmat(sqrt(sum(azidtdis.^2,2)),1,size(azidtdis,2));
    %  normalize each column to unit -- anzhao lie guiyihua
    % azidtdis = azidtdis./repmat(sqrt(sum(azidtdis.^2,1)),size(azidtdis,1),1);
    azidtdis = azidtdis';

    contourf(0:1:24,0:9:360,azidtdis,60,'Linestyle','none');
    
        %% plot
        i=1;
    for lt=0:1:24
            Hinds = azimuth1(find(localtime1>lt-0.5 ...
             & localtime1<=lt+0.5));
            Hinds_s = sind(Hinds);
            Hinds_c = cosd(Hinds);
            azi_mn_s(i) = mean(Hinds_s,"omitnan");
            azi_mn_c(i) = mean(Hinds_c,"omitnan");
            
            % calcu: alpha
            azi_mn_rad(i) = atan2(azi_mn_s(i) ,azi_mn_c(i) );
            azi_mn_deg(i) = rad2deg(azi_mn_rad(i));
            
            if azi_mn_deg(i) < 0
                azi_mn_deg(i)  = azi_mn_deg(i)  + 360;
            end
            i=i+1;
    end
    
    n_0 = find (0<azi_mn_deg & azi_mn_deg <10);
%     hold on;scatter(0:1:24,azi_mn_deg, 30, 'filled', 'MarkerFaceColor', 'k');
    hold on;plot(0:1:n_0-2,azi_mn_deg(1:1:n_0-2+1),'Color','k','LineWidth',2);
    hold on;plot(n_0-1:1:24,azi_mn_deg(n_0-1+1:1:24+1),'Color','k','LineWidth',2);

    %% label
    colormap(gca,cmap); clim([0 1]);
    ylim([0 360]);
    ytick(0:90:360)
    yminortick(0:30:360)

    xlim([0 24]);
    xtick(0:6:24)
    xminortick(0:2:24)

    %%%% LABELS
    xlabel('Local time (hour)');

end