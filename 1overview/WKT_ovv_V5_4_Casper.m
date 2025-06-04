



% WKT_ovv_V5

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
% load('D:\0lrn\0Res\Data\MCMR_MetData_2014_2023.mat');
% addpath 'D:\0lrn\0Res\Functions'
load('../../Data/MCMR_MetData_2014_2023.mat');
addpath '../../Functions'


site = 'MC'; % Mengcheng

% CAN or CAN NOT choose an equinox(zhouyepingfen) day
selected_day = datenum(2022,10,02);
inds = find(mcest>24*(selected_day- datenum(2014,01,01)) & mcest<24*(selected_day+1- datenum(2014,01,01)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTY PLOT!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; 
whitefig;
%-------------------------------------------------------

rows = 4; cols = 4;
figpos([0.35 1]);
vert_gap = 0.06*sqrt(5/rows); horz_gap = 0.14*sqrt(5/rows);
lower_marg = 0.05*sqrt(5/rows); upper_marg = 0.05*sqrt(5/rows);
left_marg = 0.064*sqrt(10/cols); right_marg = 0.018*sqrt(10/cols);

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 19; % fontsize
sz = .7*fs; % scattersize
mks = .8*fs; % markersize of meteor map

% barcolor = '#fbe3a4'; % mi yellow
barcolor = '#c0a164'; % brown
% barcolor = '#f5b265';
% barcolor = '#47b4a8';
dashcolor= '#dfdfdf';

clen = 30;
cmap = nph_saturate(cbrew('nph_RdYlBuGrey',clen),1.0);
% cmap = nph_saturate(cbrew('nph_CyclicRainbow',clen),0.8);

% axs = {...
%     [1:4, 5:8]...
%     9:10,...
%     11:12,...
%     13:14,...
%     15:16,...
%     17:20};

axs = {...
    [1:4, 5:8]...
    9:10,...
    11:12,...
    13:14,...
    15:16};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Meteors against Zenith Angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)

subplot(rows,cols,axs{5})

hold on; H = histogram(mczenith(inds),0:2.5:90,'visible','off');
hold on; bar(H.BinEdges(1:end-1)+1,H.Values,'barwidth',1,'facecolor',barcolor);
grid on;box on;

xlim([0 90])
xtick(0:30:90)
xminortick(0:5:90)

ylim([0 1100])
ytick(0:300:900)
yminortick(0:100:1100)

%%%% add zenith angle limits: % % %% why 15 to 65??
% % hold on; plot([15 15],axx.YLim,'linest','--','color',[0 0 0],'linewi',.1*fs.5);
% % hold on; plot([65 65],axx.YLim,'linest','--','color',[0 0 0],'linewi',.1*fs.5);

%%%% LABELS
ylabel('Meteor number')
xlabel('Zenith angle (degree)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Meteors against altitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)
subplot(rows,cols,axs{3})

hold on;
H=histogram(mch(inds),70:1:110,'visible','off');

hold on;
barh(H.BinEdges(1:end-1)+1,H.Values,'barwidth',1,'facecolor',barcolor);

% fit and plot fitted curve

[mcm,mcn]=hist(mch(inds),70:1:110);
[htwidth, htpeak, htA]=mygaussfit(mcn,mcm);
mcy= htA * exp( -(mcn-htpeak).^2 / (2*htwidth^2) );
hold on;plot(mcy,mcn,'color','r','linestyle','--','linewidth',.3*fs);

hpk=num2str(htpeak);hpk=hpk(1:5);
text(402,80,['\mu=' hpk ' km'], 'FontSize',fs,'fontweight','bold');
hwt=num2str(htwidth);hwt=hwt(1:5);
text(402,75,['\sigma=' hwt ' km'], 'FontSize',fs,'fontweight','bold');


grid on;box on;
xlim([0 1000])
xtick((0:300:900))
xminortick((0:100:1000))

ylim([70 110])
ytick(70:10:120)
yminortick(70:5:110)

%%%% LABELS
xlabel('Meteor number')
% xlabel('Meteor number × 10^2')
ylabel('Altitude (km)')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Meteors against Azimuth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(rows,cols,axs{4})

hold on; H1 = histogram(azimuth(inds),0:10:360,'visible','off');
% hold on; H2 = histogram(azimuth(~inrange(mczenith,[15 65])),0:10:360,'visible','off');

% hold on; bar(H1.BinEdges(1:end-1)+5,H1.Values,'barwidth',1,'facecolor',rgbtrip(0.85));
hold on; bar(H1.BinEdges(1:end-1)+5,H1.Values,'barwidth',1,'facecolor',barcolor);
% hold on; bar(H1.BinEdges(1:end-1)+5,(H1.Values-H2.Values),'barwidth',1,'facecolor',barcolor);

grid on;box on;
% axx = gca;
xlim([0 360])
xtick(0:90:360)
xminortick(0:30:360)

ylim([0 500])
ytick([0 100 300 500])
yminortick(0:50:500)

ylabel('Meteor number');
xlabel({'Azimuth (degree, clockwise from N)'})

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% Meteors against Ground Range
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % if(1)
% % 
% % subplot(rows,cols,axs{3})
% % 
% % hold on; H = histogram(mcrd,0:10:350,'visible','off');
% % hold on; bar(H.BinEdges(1:end-1),H.Values,'barwidth',1,'facecolor',barcolor);
% % % hold on; bar(H.BinEdges(1:end-1)+0.5,H.Values,'barwidth',1,'facecolor',barcolor);
% % grid on;box on;
% % 
% % end
% % 
% %         xlim([0 350])
% %         xtick(0:70:350)
% % 
% %         ylim([0 800])
% %         ytick(0:200:900)
% %         yminortick(0:50:800)
% % 
% % %%%% LABELS
% % ylabel('Meteor number')
% % xlabel('Ground range (km)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Meteors against Time of day(against local time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)

subplot(rows,cols,axs{2})
% cmap = nph_saturate(cbrew('nph_RdYlBuGrey',clen),1.0);
% hold on;
hold on; H = histogram(localtime(inds),0:24,'visible','off');
% hold on; bar(H.BinEdges(1:end-1)+0.5,H.Values,'barwidth',1,'facecolor',mcolor(6));
% hold on; bar(H.BinEdges(1:end-1)+0.5,H.Values,'barwidth',1,'facecolor',mcolor(6));
binedges = H.BinEdges(1:end-1);

% plot by time of day:
csteps = 0:(24/24):(24-(24/24));

% colour by time:
for i = 1:length(csteps)
    cinds = inrange(binedges,[csteps(i) csteps(i)+(24/24)]);
    hold on; bar(binedges(cinds)+0.5,H.Values(cinds), ...
        'barwidth',1,'facecolor',cmap(clen-i,:));
end

end

grid on;box on;
xlim([0 24])
xtick(0:6:24)
xminortick(0:3:24)

ylim([0 900])
ytick(0:300:900)
yminortick(0:100:900)

%%%% LABELS
ylabel('Meteor number');
xlabel('Local time (hour)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Meteor map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)
subplot(rows,cols,axs{1})

% get local time (24h):
localtime1 = localtime(inds);
mcrd1=mcrd(inds);
azimuth1=azimuth(inds);

% azimuth start from North clockwise????
xd1 = mcrd1.*cosd(90-azimuth1);
yd1 = mcrd1.*sind(90-azimuth1);
xd1(quadadd(xd1,yd1) > 500) = NaN;
yd1(quadadd(xd1,yd1) > 500) = NaN;

for i = 1:length(csteps)
    cinds = inrange(localtime1,[csteps(i) csteps(i)+(24/24)]);
    markerspec = {'marker','.','markersize',mks,'color',cmap(clen-i,:)};
    hold on; plot3(xd1(cinds),yd1(cinds),rand(size(xd1(cinds))),'linest','none',markerspec{:})
end

axis square
% grid on;box on;
axx = gca;
axx.XColor = 'none';
axx.YColor = 'none';

xlim(pm(350))
ylim(pm(350))
% % % COASTLINE:
% % %   if(0.0)
% % %         clat = -54.283715; clon = -36.495695;
% % %         %%%% COASTLINE:
% % %         [LatBox,~] = reckon(clat,clon,km2deg(abs(axx.YLim)),[0 180]);
% % %         [~,LonBox] = reckon(clat,clon,km2deg(abs(axx.XLim)),[90 270]);
% % %         LonBox = sort(LonBox); LatBox = sort(LatBox);
% % %         % LonBox = clon + pm(30); LatBox = clat + pm(20);
% % %         
% % %         C = nph_draw_coastline([LonBox(1) LatBox(1) ; LonBox(2) LatBox(2)],0,0,'noplot','color','k');
% % %         for q = 1:length(C)
% % %             if length(C(q).Lon) > 1
% % %                 [d,az] = distance(clat,clon,C(q).Lat,C(q).Lon);
% % %                 hold on; plot3(deg2km(d).*sind(az),deg2km(d).*cosd(az),10*ones(size(d)),'color',[1 1 1 1],'linewi',3);
% % %                 hold on; plot3(deg2km(d).*sind(az),deg2km(d).*cosd(az),10*ones(size(d)),'color',[0 0 0 1],'linewi',2);
% % %             end
%             end
%        end

%%%% let's draw some range rings and lines:
gridlinespec = {'linewi',.1*fs,'color',dashcolor};
% gridlinespec = {'linewi',.1*fs,'color',[0.15 0.15 0.15 0.15]};
hold on; plot(1.1*axx.XLim,[0 0],gridlinespec{:});
hold on; plot([0 0],1.1*axx.YLim,gridlinespec{:});

ranges = [150 250 350];
azz = 0:0.01:(2*pi);
for r = ranges
    hold on; plot3(r*sin(azz),r*cos(azz),10*ones(size(azz)),gridlinespec{:});
    % label the range lines
    textaz = 240;
    hold on; text(r*sind(textaz),r*cosd(textaz),10,[num2str(r) 'km'],'fontsize',1*fs,'horizontalalignment','center','VerticalAlignment','middle')
end

%%%% labels
hold on; text(0,1.1*axx.YLim(2),'N','fontsize',fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','bottom');
hold on; text(1.1*axx.XLim(2),0,'E','fontsize',fs,'fontweight','bold','horizontalalignment','left','VerticalAlignment','middle');
hold on; text(0,1.1*axx.YLim(1),'S','fontsize',fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
hold on; text(1.1*axx.XLim(1),0,'W','fontsize',fs,'fontweight','bold','horizontalalignment','right','VerticalAlignment','middle');

%%%% date
hold on; text(450*sind(240),450*cosd(230),[datestr((selected_day),'mmm dd, yyyy')],'fontsize',fs,'horizontalalignment','right','VerticalAlignment','middle');

%%%% number
hold on; text(450*sind(240),450*cosd(240),['# = ' num2str(length(inds))],'fontsize',fs,'horizontalalignment','right','VerticalAlignment','middle');

%%%% Mengcheng label
hold on; text(0,0,10,{[' ' upper(site)],'/'},'fontsize',fs,'fontweight','bold','horizontalalignment','left','VerticalAlignment','bottom');
set(gca,'layer','top','clipping','off')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Meteors Number against Year
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0)
subplot(rows,cols,axs{6})
bar(number./1000,'barwidth',1,'facecolor',barcolor);

%%%% LIMITS and ticks
ylim([0 24000]./1000)
ytick((0:6000:24000)./1000)
yminortick((0:1000:24000)./1000)

% day range
years = 2014:2023;

xticknyrs(gca, years, 'none');
xticknyrslabel(gca, years,fs);

%%%% LABELS
ylabel('Meteor number x 10^3');

grid on;box on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GLOBAL FORMATTING FOR ALL AXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ax = 1:5
    subplot(rows,cols,axs{ax})
    axx = gca;

    setfont(fs);
    set(gca,'linewi',.1*fs,'tickdir','out');

    switch ax
        case 1
            hold on; nph_text([0.1 0.88],['(' alphabet(ax) ')'],'fontsize',1.5*fs); %,'textborder','w');
        case {5 4}
            hold on; nph_text([0.03 0.849],['(' alphabet(ax) ')'],'fontsize',1.5*fs); %,'textborder','w');
        case {2 3}
            hold on; nph_text([0.87 0.849],['(' alphabet(ax) ')'],'fontsize',1.5*fs); %,'textborder','w');
        % case {6}
        %     hold on; nph_text([-0.015 0.849],['(' alphabet(ax) ')'],'fontsize',1.5*fs); %,'textborder','w');
    end

%     switch ax
%         case 1
%             hold on; nph_text([0.1 0.88],['(' alphabet(ax) ')'],'fontsize',1.5*fs); %,'textborder','w');
%         case {5 2 3 4}
%             hold on; nph_text([0.000 1.05],['(' alphabet(ax) ')'],'fontsize',1.5*fs); %,'textborder','w');
%         case {6}
%             hold on; nph_text([-0.035 1],['(' alphabet(ax) ')'],'fontsize',1.5*fs); %,'textborder','w');
%     end
    
end

time_sys2=clock;
timerun_21=etime(time_sys2,time_sys1);
time21=num2str(floor(timerun_21));

error
%% EXPORT FIG ==============================================================
% Set renderer to "painters"
set(gcf, 'Renderer', 'painters');
figure_name = ['Fig1'];
print(gcf,figure_name,'-djpeg','-r600');

% print(gcf,['WKT_ovv_V4_t=',time21,'s'],'-dpng','-r600');
% saveas(gcf,['WKT_ovv_V4_t=',time21,'s'],'svg');
% saveas(gcf,['Figure1'],'epsc');
disp('Figure Saved.')

function xticknyrs(gca, years, dashcolor)
daysn = datenum(years(end),12,31)-datenum(years(1),01,01)+1;
xlim([0 daysn]);

hold on;grid on; axx=gca;

% Just show the Januarys as major ticks
axx.XTick = datenum(years,01,01)-datenum(years(1),01,01);
axx.XMinorTick = 'on';
axx.XAxis.MinorTickValues = datenum(min(years),1:((length(years)+1)*12),01)-datenum(years(1),01,01);
datetick('x','m','keepticks','keeplimits');
axx.XTickLabel = {};

% add dashed lines to the Januarys except the first Jan
dashday = datenum(years(2):years(end),01,01)-datenum(years(1),01,01);
for xt = dashday
    if inrange(xt,xlim)% && ~any(xt == axx.XTick)
    hold on; plot([xt xt],axx.YLim,'linest','--','color',dashcolor,'linewi',axx.LineWidth);
%     text(xt,-500./1000,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
    end
end

end


function xticknyrslabel(gca, years, fs)
% xlim([0 datenum((max(years)+1),01,01)-datenum(min(years),01,01)]);
xlim([0 datenum((max(years)),08,01)-datenum(min(years),01,01)]);

hold on;grid on; axx=gca;

% Januarys and Junes as labels using text:
ytix1 = min(gca().YLim)-(max(gca().YLim)-min(gca().YLim))*0.05;
xtixJan = datenum(years(1),1:12:(length(years)*12+1),2)-datenum(years(1)-1,12,31);
for xt = xtixJan
    if inrange(xt,xlim) % && ~any(xt == axx.XTick)
        mn = monthname(month(xt),'mmm');
        text(xt,ytix1,mn,'fontsize',fs,'horizontalalignment','center','VerticalAlignment','top');
    end
end

xtixJul = datenum(years(1),7:12:(length(years)*12+7),2)-datenum(years(1)-1,12,31);
for xt = xtixJul
    if inrange(xt,xlim) % && ~any(xt == axx.XTick)
        mn = monthname(month(xt),'mmm');
        text(xt,ytix1,mn,'fontsize',fs,'horizontalalignment','center','VerticalAlignment','top');
    end
end

% % other months as labels using text:
% xtix = datenum(years(1),1:(length(years)*12),15)-datenum(2014,01,01);
% for xt = xtix
%     if inrange(xt,xlim) && ~any(xt == axx.XTick)
%         mn = monthname(month(xt));
%         text(xt,-500./1000,mn(1),'fontsize',fs,'horizontalalignment','center','VerticalAlignment','top');
%     end
% end

% Years as lower numbers:
ytix2 = min(gca().YLim)-(max(gca().YLim)-min(gca().YLim))*0.17;
xtix2 = datenum(years,07,01)-datenum(years(1),01,01);
for xt = xtix2
    if inrange(xt,xlim) && ~any(xt == axx.XTick)
        yrstr = datestr(xt,'yyyy');
        yrstr = max(double(yrstr))+years(1)-56+8;
        yrstr = int2str(yrstr);
%         bold number
%         text(xt,ytix2,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
        text(xt,ytix2,yrstr,'fontsize',1*fs,'horizontalalignment','center','VerticalAlignment','top');
    end
end

end