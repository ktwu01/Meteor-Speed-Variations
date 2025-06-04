



% WKTaz_sp_vec_eli_V3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ AND TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % num: start from 2014-01-01
% % % % % % est: start from 2014-01-01 0:00 UT

% % % % % % 24*(datenum(2022,12,31)- datenum(2014,01,01)+1)=max(mcest) 
% % % % % % min(mcest)=1.7575e+03 (73 days * 24 hours=1752)

% % % % % % recorded meteors: start from 74rd day 5.5 hour of 2014


% load original data
clear all;clc;close all; 
% time_sys1=clock;
% load('D:\000learning\00research\00DataAndCode\Mengcheng_wkt_Nhindley\OptionalMMRdata2014_2022.mat');
% load('D:\000learning\00research\00DataAndCode\Mengcheng_wkt_Nhindley\ids2014_2022.mat');
% addpath 'D:\000Learning\00Research\00DataAndCode\Mengcheng_wkt_Nhindley\Functions'
% load('../../Data/OptionalMMRdata2014_2022.mat');
load('../../Data/MCMR_MetData_2014_2023.mat');
load('../../Data/MCMR_ids_2022.mat');
addpath '../../Functions'

% day range
% years = 2014:2022;

fs = 20; % fontsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTY PLOT!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; 
whitefig;
% %-------------------------------------------------------
% 
rows = 2; cols = 2;
figpos([0.6 1]);

% rows = 1; cols = 4;
% figpos([1 0.6]);

vert_gap = 0.10*sqrt(5/rows); horz_gap = 0.04*sqrt(5/rows);
lower_marg = 0.025*sqrt(5/rows); upper_marg = 0.05*sqrt(5/rows);
left_marg = 0.02*sqrt(10/cols); right_marg = 0.02*sqrt(10/cols);
subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

subplot(rows,cols,1); met_vec_24(id_se,speed,localtime,azimuth,fs); % title({'Jan 4';'';'';''})
subplot(rows,cols,2); met_vec_24(id_pr,speed,localtime,azimuth,fs); % title({'Apr 4';'';'';''})
subplot(rows,cols,3); met_vec_24(id_ae,speed,localtime,azimuth,fs); % title({'Jul 4';'';'';''})
subplot(rows,cols,4); met_vec_24(id_ap,speed,localtime,azimuth,fs); % title({'Oct 2';'';'';''})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GLOBAL FORMATTING FOR ALL AXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MonthDay = {' Jan 4';' Apr 4';' Jul 4';' Oct 2'};

for ax = 1:rows*cols
    subplot(rows,cols,ax);
    axx = gca;
    setfont(fs);
    % set(gca,'linewi',1,'tickdir','out');

%     switch ax
%         case 1
            % hold on; nph_text([0.015 0.88],strcat('(', alphabet(ax), ')', MonthDay{ax}),'fontsize',1.5*fs); %,'textborder','w');
            % hold on; nph_text([0.015 0.93],strcat('(', alphabet(ax), ')', MonthDay{ax}),'fontsize',1.5*fs); %,'textborder','w');
            hold on; nph_text([0.015 0.98],strcat('(', alphabet(ax), ')', MonthDay{ax}),'fontsize',1.5*fs); %,'textborder','w');
%         case {5 2 4}
%             hold on; nph_text([0.05 0.87],['(' alphabet(ax) ')'],'fontsize',1.5*fs); %,'textborder','w');
%         case {3}
%             hold on; nph_text([0.85 0.87],['(' alphabet(ax) ')'],'fontsize',1.5*fs); %,'textborder','w');
%         case {6}
%             hold on; nph_text([-0.015 0.87],['(' alphabet(ax) ')'],'fontsize',1.5*fs); %,'textborder','w');
%     end

end

% time_sys2=clock;
% timerun_21=etime(time_sys2,time_sys1);
% time21=num2str(floor(timerun_21));


error

%% EXPORT FIG ==============================================================
figure_name = 'Fig7B';
print(gcf,figure_name,'-djpeg','-r600');

% saveas(gcf,['WKTaz_sp_vec_eli_V2_t=',time21,'s'],'png')
% saveas(gcf,['WKTaz_sp_vec_eli_V2_t=',time21,'s'],'tif')
% saveas(gcf,['WKTaz_sp_vec_eli_V2_t=',time21,'s'],'fig')
disp('Figure Saved.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vec: meteor as a function of azimuth and speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% E.G. inds0 = id_pr;


function met_vec_24(inds0,speed,localtime,azimuth,fs)

% fs = input; % fontsize
dashcolor= '#dfdfdf';

axis square
% grid on;box on;
axx = gca;
axx.XColor = 'none';
axx.YColor = 'none';


% limit matters ! it depends on the mode value of the max vector
lim0 = 15; % km/s

xlim(pm(lim0))
ylim(pm(lim0))

%%%% let's draw some range rings and lines:
gridlinespec = {'linewi',1,'color',dashcolor};
% gridlinespec = {'linewi',1,'color',[0.15 0.15 0.15 0.15]};
hold on; plot(1.1*axx.XLim,[0 0],gridlinespec{:});
hold on; plot([0 0],1.1*axx.YLim,gridlinespec{:});

ranges = lim0/3: lim0/3: lim0;
azz = 0:0.01:(2*pi);
for r = ranges
    hold on; plot3(r*sin(azz),r*cos(azz),10*ones(size(azz)),gridlinespec{:});
    % label the range lines
    textaz = 240;
    hold on; text(r*sind(textaz),r*cosd(textaz),10,[num2str(r) ' km/s'],'fontsize',1*fs,'horizontalalignment','center','VerticalAlignment','middle')
end

%%%% labels
hold on; text(0,1.1*axx.YLim(2),{'N';''},'fontsize',fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','bottom');
hold on; text(1.1*axx.XLim(2),0,'      E','fontsize',fs,'fontweight','bold','horizontalalignment','left','VerticalAlignment','middle');
hold on; text(0,1.1*axx.YLim(1),{'S'},'fontsize',fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
hold on; text(1.1*axx.XLim(1),0,'W      ','fontsize',fs,'fontweight','bold','horizontalalignment','right','VerticalAlignment','middle');

%%%% date
% hold on; text(450*sind(240),450*cosd(230),['2022,10,02'],'fontsize',fs,'horizontalalignment','right','VerticalAlignment','middle');

%%%% Mengcheng label
site = 'MC'; % Mengcheng
hold on; text(0,0,10,{[' ' upper(site)],'/'},'fontsize',fs,'fontweight','bold','horizontalalignment','left','VerticalAlignment','bottom');
% set(gca,'layer','top','clipping','off')

metvec_az_sp(inds0,	0	,	1	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	1	,	2	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	2	,	3	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	3	,	4	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	4	,	5	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	5	,	6	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	6	,	7	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	7	,	8	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	8	,	9	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	9	,	10	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	10	,	11	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	11	,	12	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	12	,	13	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	13	,	14	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	14	,	15	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	15	,	16	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	16	,	17	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	17	,	18	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	18	,	19	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	19	,	20	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	20	,	21	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	21	,	22	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	22	,	23	,speed,localtime,azimuth,fs)	
metvec_az_sp(inds0,	23	,	24	,speed,localtime,azimuth,fs)	
end


function metvec_az_sp(inds0,tmin,tmax,speed,localtime,azimuth,fs)

lct = localtime(inds0);
inds = inds0(lct > tmin & lct < tmax);
spd = speed(inds)./1000;
azimuth1=azimuth(inds);

% fs = input; % fontsize
% sz = 7; % scattersize
% mks = 6; % markersize of meteor map

% % cmap
% % clen matters! this means the max limit of speed considered
% clen = 80; cmap = flip(nph_saturate(cbrew('nph_alt_jet',clen),1.0));


%%%% let's calculate the vector!
U = mean (spd.*cosd(90-azimuth1));
V = mean (spd.*sind(90-azimuth1));
% W = sqrt(U*U+V*V);
% phi = asin(V/W);
% r = W * exp(i * phi);
hold on;%compass(r);
asf = 1;
quiver (0,0,U,V,'AutoScaleFactor',asf,'MaxHeadSize',0.5,'LineWidth',4);

set(gca,'layer','top','clipping','off')

%%%% lct
hold on; text(U*asf,V*asf,[num2str(tmin) '-' num2str(tmax)],'fontsize',fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','middle');

% [num2str(tmin) ' < LT < ' num2str(tmax)]
% % title({[num2str(tmin) ' < Local time < ' num2str(tmax)];'';''})

end


% test
% inds0 = id_ae;
% tmin = 6
% tmax = 9
% 
% lct = localtime(inds0);
% inds = inds0(lct > tmin & lct < tmax);
% spd = speed(inds)./1000;
% azimuth1=azimuth(inds);
% 
% U = mean (spd.*cosd(90-azimuth1))*50
% V = mean (spd.*sind(90-azimuth1))*50
