


%%% zeni_lt_surf_con_1_4p_V2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ AND TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % num: start from 2014-01-01
% % % % % % est: start from 2014-01-01 0:00 UT

% % % % % % 24*(datenum(2022,12,31)- datenum(2014,01,01)+1)=max(mcest) 
% % % % % % min(mcest)=1.7575e+03 (73 days * 24 hours=1752)

% % % % % % recorded meteors: start from 74rd day 5.5 hour of 2014

clear all; clc; close all; 
% time_sys1=clock;
% load original data
% load('D:\0lrn\00Res\Data\OptionalMMRdata2014_2022.mat', ...
%     'mczenith','localtime','mcest');??
% load(['D:\0lrn\00Res\Data\ids2014_2022.mat']);??
load('../../Data/MCMR_MetData_2014_2023.mat','mczenith','localtime','mcest');
load('../../Data/MCMR_ids_2022.mat');
addpath '../../Functions'

% addpath 'D:\0lrn\00Res\Functions'

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
% fs = 10; % fontsize
fs = 20; % fontsize
% sz = .6*fs; % scattersize
% mks = .3*fs; % markersize of meteor map

axs = {
    [1:8 9:16],
    17:18,
    19:20,
    21:22,
    23:24
    };

subplot(rows,cols,axs{1});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calcu surf
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
azimuth1 = mczenith(inds);
% speed1=speed(inds)/1000;

azidtdis = zeros(72+1,31);
n=1;
for lt=0:1:72
    m=1;
    for az=0:3:90
        azidtdis(n,m)=length(find(azimuth1>az-1.5 & ...
            azimuth1<az+1.5 ...
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
title(datestr((datenum(2022,10,02)),'mmm dd, yyyy'));
[xx,yy] = meshgrid (0:3:90,0:1:72);

%% 
surf(xx,yy,azidtdis,'FaceAlpha',1,'LineStyle','none');

hold on;
plot3(0:3:90,repmat(0,1,31),azidtdis(1,:),LineWidth=2,Color='r');
hold on;
plot3(repmat(0,1,72+1),0:1:72,azidtdis(:,1),LineWidth=2,Color='r');

xlim([0 90.0]);
xtick(0:15:90.0)
xminortick(0:3:90)

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
xlabel('Zenith (degree)');

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

    axx = gca;
    setfont(fs);
    set(gca,'linewi',0.5,'tickdir','out');
    hold on; nph_text([-0.015 0.87],['(' alphabet(1) ')'],'fontsize',1.5*fs);

subplot(rows,cols,axs{2});azidtdis1=az_lt_contourf(id_pr,localtime, mczenith, 1);title('Jan 4')
    ylabel('Zenith (degree)');
    axx = gca;
    setfont(fs);
    set(gca,'linewi',0.5,'tickdir','out');
    hold on; nph_text([-0.17 1.1],['(' alphabet(2) ')'],'fontsize',1.5*fs); %,'textborder','w');

subplot(rows,cols,axs{3});azidtdis2=az_lt_contourf(id_se,localtime, mczenith, 2);title('Apr 4')
    axx = gca;
    setfont(fs);
    set(gca,'linewi',0.5,'tickdir','out');
    hold on; nph_text([-0.17 1.1],['(' alphabet(3) ')'],'fontsize',1.5*fs); %,'textborder','w');

subplot(rows,cols,axs{4});azidtdis3=az_lt_contourf(id_ap,localtime, mczenith, 3);title('Jul 4')
    axx = gca;
    setfont(fs);
    set(gca,'linewi',0.5,'tickdir','out');
    hold on; nph_text([-0.17 1.1],['(' alphabet(4) ')'],'fontsize',1.5*fs); %,'textborder','w');

subplot(rows,cols,axs{5});azidtdis4=az_lt_contourf(id_ae,localtime, mczenith, 4);title('Oct 2')
    cbar= nph_colorbar;
    axx = gca;
    setfont(fs);
    set(gca,'linewi',0.5,'tickdir','out');
    hold on; nph_text([-0.17 1.1],['(' alphabet(5) ')'],'fontsize',1.5*fs); %,'textborder','w');

% for ax = 1:5
%     subplot(rows,cols,axs{ax});
%     axx = gca;
%     setfont(fs);
%     set(gca,'linewi',0.5,'tickdir','out');
%     switch ax
%         case 1
%         hold on; nph_text([-0.015 0.87],['(' alphabet(ax) ')'],'fontsize',1.5*fs);
%         case {2 3 4 5}
%         hold on; nph_text([-0.17 1.1],['(' alphabet(ax) ')'],'fontsize',1.5*fs); %,'textborder','w');
%     end
% end

% time_sys2=clock;
% timerun_21=etime(time_sys2,time_sys1);

error

%% EXPORT FIG ==============================================================

figure_name = ['Fig5'];
print(gcf,figure_name,'-djpeg','-r600');

% saveas(gcf,'WKTzeni_lt_surf_con_1_4p_V1','png')
% saveas(gcf,'WKTzeni_lt_surf_con_1_4p_V1','epsc')
% saveas(gcf,'WKTzeni_lt_surf_con_1_4p_V1','fig')
disp('Fig:WKTzeni_lt_surf_con_1_4p_V1 Saved.')

function azidtdis = az_lt_contourf(id_ae,localtime, azimuth, px)
    % colormap
    clen=100;cmap = nph_saturate(cbrew('nph_alt_jet',clen),1.0);

    % calcu contourf
    localtime1=localtime(id_ae);
    azimuth1 = azimuth(id_ae);
    
    azidtdis = zeros(25,31);
    n=1;
    
    for lt=0:1:24
        m=1;
        for az=0:3:90
            azidtdis(n,m)=length(find(azimuth1>az-1.5 & ...
                azimuth1<az+1.5 ...
                & localtime1>lt-0.5 & localtime1<lt+0.5));
            m=m+1;
        end
        n=n+1;
    end
    
    % %  normalize each row to unit
    for lt = (0:1:24)+1
    azidtdis(lt,:) = azidtdis(lt,:)/max(azidtdis(lt,:));
    end

    %% calculate the peak
        i=1;
    for lt=0:1:24
            Hinds = azimuth1(find(localtime1>lt-0.5 ...
             & localtime1<=lt+0.5));
            %% FIT start
            H = histogram (Hinds,0:1:90,'Visible','off');
            x = H.BinEdges(1:end-1)+.5;
            y = H.Values;
            y = y/max(y);

            f = fit(x.',y.','gauss2', ...% a1, mu1, sig1, a2, mu2, sig2
            'Lower',[0.5, 35.0, 15,-1, 60, 7], ...
            'Start',[1.0, 60.0, 20,-0.5, 75, 12], ... 
            'Upper',[1.5, 70.0, 30, 0,   85, 15]);
            yi= f.a1*exp(-((x-f.b1)/f.c1).^2) + f.a2*exp(-((x-f.b2)/f.c2).^2);
            [~, idx] = max(yi);
            mu = x(idx);

            % Find peak value
            peak_value(i) = mu;
            i=i+1;
    end

    azidtdis = azidtdis';

    hold on; contourf(0:1:24,0:3:90,azidtdis,100,'Linestyle','none');
    hold on; plot(0:1:24,peak_value,'Color','k','LineWidth',2);
    colormap(gca,cmap); clim([0 1]);

    %% label
    ylim([0 90]);
    ytick(0:15:90)
    yminortick(0:3:90)

    xlim([0 24]);
    xtick(0:6:24)
    xminortick(0:2:24)

    %%%% LABELS
    xlabel('Local time (hour)');

end
