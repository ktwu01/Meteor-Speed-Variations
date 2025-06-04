



% WKTsp_az_map_V2

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
load('../../Data/OptionalMMRdata2014_2022.mat','ndays'); % ndays = 3214
load('../../Data/MCMR_MetData_2014_2023.mat');
load(['../../Data/MCMR_ids_2022.mat']);
addpath '../../Functions'

% day range
% years = 2014:2022;
fs = 9; % fontsize
% fs = 18; % fontsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTY PLOT!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all;
whitefig;
%-------------------------------------------------------

rows = 4; cols = 4;
figpos([0.7 1]);
vert_gap = 0.05*sqrt(5/rows); horz_gap = 0.05*sqrt(5/rows);
lower_marg = 0.05*sqrt(5/rows); upper_marg = 0.055*sqrt(5/rows);
left_marg = 0.05*sqrt(10/cols); right_marg = 0.03*sqrt(10/cols);

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

% color scale by speeD:
csteps = 0:(80/80):(80-(80/80));

Create_Q = 1

Create_R = 1

subplot(rows,cols,1)
metmap_az_sp(id_ae,0,3,speed,localtime,azimuth,mcrd,csteps,fs,Create_Q,Create_R)

subplot(rows,cols,2)
metmap_az_sp(id_ae,3,6,speed,localtime,azimuth,mcrd,csteps,fs,Create_Q,Create_R)

subplot(rows,cols,3)
metmap_az_sp(id_ae,6,9,speed,localtime,azimuth,mcrd,csteps,fs,Create_Q,Create_R)

subplot(rows,cols,4)
metmap_az_sp(id_ae,9,12,speed,localtime,azimuth,mcrd,csteps,fs,Create_Q,Create_R)

subplot(rows,cols,5)
met_sp_hist(id_ae,0,3,speed,localtime,csteps,ndays,fs)
% %%%% LABELS
ylabel('Meteor Number');

subplot(rows,cols,6)
met_sp_hist(id_ae,3,6,speed,localtime,csteps,ndays,fs)
subplot(rows,cols,7)
met_sp_hist(id_ae,6,9,speed,localtime,csteps,ndays,fs)
subplot(rows,cols,8)
met_sp_hist(id_ae,9,12,speed,localtime,csteps,ndays,fs)


subplot(rows,cols,9)
metmap_az_sp(id_ae,12,15,speed,localtime,azimuth,mcrd,csteps,fs,Create_Q,Create_R)

subplot(rows,cols,10)
metmap_az_sp(id_ae,15,18,speed,localtime,azimuth,mcrd,csteps,fs,Create_Q,Create_R)

subplot(rows,cols,11)
metmap_az_sp(id_ae,18,21,speed,localtime,azimuth,mcrd,csteps,fs,Create_Q,Create_R)

subplot(rows,cols,12)
metmap_az_sp(id_ae,21,24,speed,localtime,azimuth,mcrd,csteps,fs,Create_Q,Create_R)

subplot(rows,cols,13);
met_sp_hist(id_ae,12,15,speed,localtime,csteps,ndays,fs)
% %%%% LABELS
ylabel('Meteor Number');
xlabel('Speed (km/s)')

subplot(rows,cols,14)
met_sp_hist(id_ae,15,18,speed,localtime,csteps,ndays,fs)
xlabel('Speed (km/s)')

subplot(rows,cols,15)
met_sp_hist(id_ae,18,21,speed,localtime,csteps,ndays,fs)
xlabel('Speed (km/s)')

subplot(rows,cols,16)
met_sp_hist(id_ae,21,24,speed,localtime,csteps,ndays,fs)
xlabel('Speed (km/s)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GLOBAL FORMATTING FOR ALL AXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ax = 1:rows*cols
    subplot(rows,cols,ax);
    axx = gca;
    setfont(fs);
    set(gca,'linewi',1,'tickdir','out');

    %     switch ax
    %         case 1
    hold on; nph_text([0.021 0.85],['(' alphabet(ax) ')'],'fontsize',1.5*fs); %,'textborder','w');
    %         case {5 2 4}
    %             hold on; nph_text([0.05 0.87],['(' alphabet(ax) ')'],'fontsize',1.5*fs); %,'textborder','w');
    %         case {3}
    %             hold on; nph_text([0.85 0.87],['(' alphabet(ax) ')'],'fontsize',1.5*fs); %,'textborder','w');
    %         case {6}
    %             hold on; nph_text([-0.015 0.87],['(' alphabet(ax) ')'],'fontsize',1.5*fs); %,'textborder','w');
    %     end

end

time_sys2=clock;
timerun_21=etime(time_sys2,time_sys1);
time21=num2str(floor(timerun_21));

error
%% EXPORT FIG ==============================================================
% Set renderer to "painters"
set(gcf,'Renderer', 'painters');
figure_name = ['Fig7'];
print(gcf,figure_name,'-djpeg','-r600');
% print(gcf,figure_name,'-dpng','-r600');
% saveas(gcf,['Figure7'],'epsc');
% saveas(gcf,figure_name,'svg');
disp('Figure Saved.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Meteors against speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function met_sp_hist(inds0,tmin,tmax,speed,localtime,csteps,ndays,fs)
% fs = this function's input; % fontsize
% dashcolor= '#dfdfdf';

% cmap
% clen matters! this means the max limit of speed considered
clen = 80; cmap = flip(nph_saturate(cbrew('nph_alt_jet',clen),1.0));

lct = localtime(inds0);
spd = speed(inds0)./1000;
inds = inds0(lct > tmin & lct < tmax & spd <80);
spd = speed(inds)./1000;

H = histogram(spd,0:80,'visible','off');
binedges = H.BinEdges(1:end-1);

% colour by time:
for i = 1:length(csteps)-1
    cinds = inrange(binedges,[csteps(i) csteps(i)+(80/80)]);
    hold on; bar(binedges(cinds)+0.5,H.Values(cinds)./ndays, ...
        'barwidth',1,'facecolor',cmap(clen-i,:));
end

grid on;box on;
xlim([0 80])
xtick(0:20:80)
xminortick(0:10:80)

% ylim([0 900])
% ytick(0:300:900)
% yminortick(0:100:900)

% %%%% LABELS
% ylabel('Meteor Number');
% xlabel('Speed (km/s)')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Map: meteor as a function of azimuth and speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% E.G. inds0 = id_pr;

function metmap_az_sp(inds0,tmin,tmax,speed,localtime,azimuth,mcrd,csteps,fs,Create_Q,Create_R)

% fs = input; % fontsize
% sz = fs*7/9; % scattersize
mks = fs*.5; % markersize of meteor map
% sz = 7; % scattersize
% mks = 6; % markersize of meteor map
dashcolor= '#dfdfdf';

% cmap
% clen matters! this means the max limit of speed considered
clen = 80; cmap = flip(nph_saturate(cbrew('nph_alt_jet',clen),1.0));

lct = localtime(inds0);
inds = inds0(lct > tmin & lct < tmax);
spd = speed(inds)./1000;

mcrd1=mcrd(inds);
azmu=azimuth(inds);

% azimuth start from North clockwise????
xd1 = mcrd1.*cosd(90-azmu);
yd1 = mcrd1.*sind(90-azmu);
xd1(quadadd(xd1,yd1) > 500) = NaN;
yd1(quadadd(xd1,yd1) > 500) = NaN;

for i = 1:length(csteps)-1
    cinds = inrange(spd,[csteps(i) csteps(i)+(80/80)]);
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
gridlinespec = {'linewi',1,'color',dashcolor};
% gridlinespec = {'linewi',1,'color',[0.15 0.15 0.15 0.15]};
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
set(gca,'layer','bottom','clipping','off')

%%%% labels
hold on; text(0,1.1*axx.YLim(2),'N','fontsize',fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','bottom');
hold on; text(1.1*axx.XLim(2),0,'E','fontsize',fs,'fontweight','bold','horizontalalignment','left','VerticalAlignment','middle');
hold on; text(0,1.1*axx.YLim(1),'S','fontsize',fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
hold on; text(1.1*axx.XLim(1),0,'W','fontsize',fs,'fontweight','bold','horizontalalignment','right','VerticalAlignment','middle');

%%%% date
% hold on; text(450*sind(240),450*cosd(230),['2022,10,02'],'fontsize',fs,'horizontalalignment','right','VerticalAlignment','middle');

%%%% angle: the degree means clocwise from north
%%%% lct
textaz1 = 240; % degree
hold on; text(450*sind(textaz1),450*cosd(textaz1),[num2str(tmin) ' < LT < ' num2str(tmax)],'fontsize',fs,'fontweight','bold','horizontalalignment','right','VerticalAlignment','middle');

%%%% number
textaz2 = 230; % degree
hold on; text(490*sind(textaz2),460*cosd(textaz2),['# = ' num2str(length(inds))],'fontsize',fs,'fontweight','bold','horizontalalignment','right','VerticalAlignment','middle');
set(gca,'layer','bottom','clipping','off')


%%%% Mengcheng label
site = 'MC'; % Mengcheng
hold on; text(0,0,10,{[' ' upper(site)],'/'},'fontsize',fs,'fontweight','bold','horizontalalignment','left','VerticalAlignment','bottom');
set(gca,'layer','bottom','clipping','off')

if (Create_Q)
    %%%% let's calculate the Q vector!
    Color_Q = "#D95319"; % color of Q and its scale bar
    sf_Q = 53; % Scale factor
    Q_U = mean (spd.*cosd(90-azmu))*sf_Q;
    Q_V = mean (spd.*sind(90-azmu))*sf_Q;
    Q = sqrt(Q_U*Q_U+Q_V*Q_V)/sf_Q;
    % hold on;
    % compass(r);
    quiver (0,0,Q_U,Q_V,'AutoScale','off','MaxHeadSize',0.5,'LineWidth',4,'Color',Color_Q);

    %%%% Add a scale bar for Q, on top of the dots
    % hold on; % compass(r);
    scale_bar_x_Q = 350;
    exScale_Q = 5; % exScale_Q km/s
    bar_height_Q = sf_Q * exScale_Q; % km in the map
    quiver (scale_bar_x_Q,0,0,bar_height_Q,'AutoScale','off','MaxHeadSize',0.5,'LineWidth',2,'Color',Color_Q);
    text(scale_bar_x_Q, bar_height_Q, sprintf('%d km/s', exScale_Q), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize',fs, 'Rotation', 90,'Color',Color_Q);

    % %%%% Notes to Q
    % hold on; text(U,V,['Q = ' num2str(Q)],'fontsize',fs,'fontweight','bold','horizontalalignment','left','VerticalAlignment','middle');
    % title({[num2str(tmin) ' < Local time < ' num2str(tmax)];'';''})

end


if (Create_R)

    %%%% let's calculate the R vector!
    %%%% the velocity averaging at each azimuth (named R) as a function of azimuth
    %%%% need to calculate the averaged speed at each azimuth
    %%%% and then plot these values on the map
    %%%% R do not need to have the same scale bar as Q
    sf_R = 13;
    Color_R = 'k'; % color of R and its scale bar
    % Calculate and plot R vectors
    azimuth_bins = 0:10:360; % Define azimuth bins
    R_U = zeros(size(azimuth_bins));
    R_V = zeros(size(azimuth_bins));

    for i = 1:length(azimuth_bins)-1
        bin_inds = azmu >= azimuth_bins(i) & azmu < azimuth_bins(i+1);
        if any(bin_inds)
            R_U(i) = mean(spd(bin_inds).*cosd(90-azmu(bin_inds)))*sf_R;
            R_V(i) = mean(spd(bin_inds).*sind(90-azmu(bin_inds)))*sf_R;
        end
    end

    R_U(end) = R_U(1);
    R_V(end) = R_V(1);

    % Plot R vectors end-to-end like a radar map

    for i = 1:length(azimuth_bins)-1
        bin_inds = azmu >= azimuth_bins(i) & azmu < azimuth_bins(i+1);
        if any(bin_inds)
            x_start = R_U(i);
            y_start = R_V(i);
            x_end = R_U(i+1);
            y_end = R_V(i+1);
            hold on; plot([x_start x_end], [y_start y_end], 'LineWidth', 2, 'Color', Color_R);
        end
    end

    % % % Calculate magnitude of R vectors
    % % R_magnitude = sqrt(R_U.^2 + R_V.^2);
    % % % Visualize the histogram of the magnitudes
    % % figure;
    % % plot(R_magnitude);
    % % title('Histogram of R Magnitudes');
    % % xlabel('Magnitude of R');
    % % ylabel('Frequency');

% R = sqrt(R_U.*R_U+R_V.*R_V)/sf_R; % calculate each R
%%%% Add a scale bar for R, on top of the dots
% hold on; % compass(r);
scale_bar_x_R = 350;
exScale_R = 30; % exScale_R km/s
bar_height_R = -sf_R * exScale_R; % km in the map
quiver (scale_bar_x_R,0,0,bar_height_R,'AutoScale','off','MaxHeadSize',0.5,'LineWidth',2,'Color',Color_R);
text(scale_bar_x_R, bar_height_R, sprintf('%d km/s', exScale_R), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize',fs, 'Rotation', 90,'Color',Color_R);

% %%%% Notes to R

end
end

% test
% inds0 = id_ae;
% tmin = 6
% tmax = 9
%
% lct = localtime(inds0);
% inds = inds0(lct > tmin & lct < tmax);
% spd = speed(inds)./1000;
% azmu=azimuth(inds);
%
% U = mean (spd.*cosd(90-azmu))*50
% V = mean (spd.*sind(90-azmu))*50
