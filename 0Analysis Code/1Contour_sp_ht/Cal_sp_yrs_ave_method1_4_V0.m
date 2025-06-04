



% Cal_sp_yrs_ave_method1_4_V0.m
% load original data
clear all;clc;close all;
addpath 'D:\0lrn\00Res\Functions'
load('D:\0lrn\00Res\Data\McMRdata_2014_2023.mat','doy','speed','number','mcest');
speed = speed/1000;
DAY_PER = datenum(2023,12,31)- datenum(2014,1,1)+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Speed contourf of each day of N years
% time use: about 7 mins*DAY_PER/108
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)

% mcspdisNyrs = zeros(41,DAY_PER+1);

% create bin edges for "date" and "speed"
dateBinEdges = 0:24:24*DAY_PER;
speedBinEdges = 0:2:80+2; %adding one more bin in the end

% Use histcounts2 to calculate bins and counts
mcspdisNyrs = histcounts2(mcest,speed,dateBinEdges,speedBinEdges);
mcspdisNyrs = mcspdisNyrs';

disp('Caculated: speed contourf data of each day of N years');

% save file
save(['D:\0lrn\00Res\Data\Mc_sp_dis_daily_2014_2023.mat'], ...
    'mcspdisNyrs');
disp('Saved: speed contourf data of each day of N years');
figure;contourf(mcspdisNyrs,100,'Linestyle','none');colorbar; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% speed contourf of an ave year
% time use: about 20 mins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)

edges_dt = 0.5:366.5; % day of year
edges_sp = 0:2:80+2; %adding one more bin in the end

% first vector is doy, second is speed
[counts,~,~] = histcounts2(doy, speed, edges_dt, edges_sp);
% mcspdisave = counts
mcspdisave = counts.'; % Transpose to match original array shape

disp('Calcued: speed contourf data of ave year');
% save file
save(['D:\0lrn\00Res\Data\Mc_sp_dis_ave_2014_2023.mat'], ...
    'mcspdisave');

figure;contourf(mcspdisave,100,'Linestyle','none');colorbar; 
disp('Saved: speed contourf data of ave year');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% speed contourf of ave year: method 4
% NOTE: did not subtract the BG
% time use: about 20 mins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)

load(['D:\0lrn\00Res\Data\Mc_sp_dis_ave_2014_2023.mat'], ...
    'mcspdisave');

figure;
clen=100; cmap = nph_saturate(cbrew('nph_Rainbow',clen),1.0);
hold on; colormap(gca,cmap);
contourf(0:365,0:2:80,mcspdisave,100,'Linestyle','none');colorbar; 
ylim([10 70])


hold on; axx=gca;
axx.XTick = datenum(2014,1:12,1)-datenum(2014,01,01);
datetick('x','m','keepticks','keeplimits')

axx.XTickLabel = {};
fs = 12

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

%% smooth period
SmtPrd = 30

for i = 1:41
    smave(i,:) = smoothdata (mcspdisave(i,:),'movmedian',SmtPrd);
end

mcspdisave = (mcspdisave-smave)./smave;


figure;
clen=100; cmap = nph_saturate(cbrew('nph_RainbowWhite',clen),0.7);
hold on; colormap(gca,cmap);
contourf(0:365,0:2:80,mcspdisave,100,'Linestyle','none');colorbar; 
ylim([20 70])

% % % % % % % %% all finally save
% % % % % % % save('D:\0lrn\00Res\Data\MMRalldis_2014_2022.mat');
% % % % % % % disp('All successfully saved!');

end