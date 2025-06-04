



%%% WKTht_sp_lt_hist_4col_V1

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
% load('../../Data/OptionalMMRdata2014_2022.mat');
% load(['../../Data/ids2014_2022.mat']);
load('../../Data/MCMR_MetData_2014_2023.mat');
load('../../Data/MCMR_ids_2022.mat');


addpath '../../Functions'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure A!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)
figure; hold all; 
whitefig;

%--------------------------------------------------------
rows = 2; cols = 4;

figpos([1 0.32*rows]);
vert_gap = 0.05*sqrt(5/rows)*1;    horz_gap = 0.026*sqrt(10/cols);
lower_marg = 0.06*sqrt(5/rows);     upper_marg = 0.02*sqrt(5/rows)*1.8;
left_marg = 0.03*sqrt(10/cols);     right_marg = 0.04*sqrt(10/cols);
subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);
%--------------------------------------------------------
% fs = 12; % fontsize
fs = 20; % fontsize

%% barcolor bar color
% Example: barcolor = '#0065a1';
dashcolor = '#a5a5a5';
barcolor = '#009BCA';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SubPlots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% id_se : inds for spring equinox March 20, 2022, 11:33 PM, Mar 13-27
% id_ss : inds for summer solstice June 21, 2022, 5:13 PM, Jun 14-28
% id_ae : inds for autumn equinox September 23, 2022, 9:04 AM, Sep 16-30
% id_ws : inds for winter solstice December 22, 2022, 5:47 AM, Dec 15-29
% id_pr : inds for Perihelion Jan-4-2022 R_earthandsun= 147105052 km
% id_ap : inds for Aphelion July-4-2022 R_earthandsun= 152098455 km

%% plot
subplot(rows,cols, 1);numdis_se1 = mch_localtime_contourf(id_pr, localtime, mch, fs);title('Jan 4');ylabel('Height (km)');
subplot(rows,cols, 2);numdis_se2 = mch_localtime_contourf(id_se, localtime, mch, fs);title('Apr 4');
subplot(rows,cols, 3);numdis_se3 = mch_localtime_contourf(id_ap, localtime, mch, fs);title('Jul 4');
subplot(rows,cols, 4);numdis_se4 = mch_localtime_contourf(id_ae, localtime, mch, fs);title('Oct 2');

% %%%% COLORBAR
% hold on;
% cbar = nph_colorbar;clim([0 1]);
% cbar.Label.String = ['N/N_{max}']; %cbar.Label.String = ['log_{10}(N/N_{max})'];
% cbar.TickDirection = 'out';
% cbar.Position = cbar.Position .* [1 1 0.8 1];
% cbar.Position = cbar.Position + [0.5*cbar.Position(3) 0 0 0];
% %             cbar.Label.String = ['\bf{' clabs{ax} '}'];
% cbar.Label.Rotation = 90;

%% plot
subplot(rows,cols, 5);numdis_se5 = speed_localtime_contourf(id_pr, localtime, speed, fs);ylabel('Speed (km/s)');
subplot(rows,cols, 6);numdis_se6 = speed_localtime_contourf(id_se, localtime, speed, fs);
subplot(rows,cols, 7);numdis_se7 = speed_localtime_contourf(id_ap, localtime, speed, fs);
subplot(rows,cols, 8);numdis_se8 = speed_localtime_contourf(id_ae, localtime, speed, fs);

% %%%% COLORBAR
% hold on;
% cbar = nph_colorbar;clim([0 1]);
% cbar.Label.String = ['N/N_{max}']; %cbar.Label.String = ['log_{10}(N/N_{max})'];
% cbar.TickDirection = 'out';
% cbar.Position = cbar.Position .* [1 1 0.8 1];
% cbar.Position = cbar.Position + [0.5*cbar.Position(3) 0 0 0];
% %             cbar.Label.String = ['\bf{' clabs{ax} '}'];
% cbar.Label.Rotation = 90;

%%%% COLORBAR
hold on;
cbar = nph_colorbar;clim([0 1]);
cbar.Label.String = ['N/N_{max}']; %cbar.Label.String = ['log_{10}(N/N_{max})'];
cbar.TickDirection = 'out';
cbar.Position = cbar.Position .* [1 1 0.8 1];
cbar.Position = cbar.Position + [0.5*cbar.Position(3) 21*cbar.Position(3) 0 0];
%             cbar.Label.String = ['\bf{' clabs{ax} '}'];
cbar.Label.Rotation = 90;


%% plot
% subplot(rows,cols, 5);mcount_localtime_hist(id_pr, localtime);ylabel('Meteor Number × 10^2');
% subplot(rows,cols, 6);mcount_localtime_hist(id_se, localtime);
% subplot(rows,cols, 7);mcount_localtime_hist(id_ap, localtime);
% subplot(rows,cols, 8);mcount_localtime_hist(id_ae, localtime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GLOBAL FORMATTING FOR ALL AXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ax = 1:rows*cols
%     subplot(rows,cols,axs{ax});
    subplot(rows,cols,ax);
    
    axx = gca;
    setfont(fs);
    set(gca,'LineWidth',.08*fs,'tickdir','out');

hold on; nph_text([0.015 0.895],['(' alphabet(ax) ')'],'fontsize',1.5*fs); 
end

% time_sys2=clock;
% timerun_21=etime(time_sys2,time_sys1);
% time21=num2str(floor(timerun_21));

error

%% EXPORT FIG ==============================================================
% Set renderer to "painters"
set(gcf, 'Renderer', 'painters');

figure_name = ['Fig6'];
print(gcf,figure_name,'-djpeg','-r600');
% print(gcf,['WKTht_sp_lt_hist_4col_V1_A_t='],'-dpng','-r600');
% saveas(gcf,['Figure6'],'epsc');
% saveas(gcf,['WKTht_sp_lt_hist_4col_V1_A_t=',time21,'s'],'svg');
disp('Figure A Saved.');
end

%% speed_localtime_contourf
function numdis_se = speed_localtime_contourf(id_se, localtime, speed, fs)
localtime_se= localtime(id_se);% get Local time (24h)
speed_se=speed(id_se)./1000;

numdis_se = zeros(41,25);
n=1;
for sp=0:2:80
    m=1;
    for hr=0:1:24
        numdis_se(n,m)=length(find(localtime_se>hr-0.5 & ...
            localtime_se<=hr+0.5 & speed_se>sp-1 & speed_se<sp+1));
        m=m+1;
    end
    n=n+1;
end

% %  normalize each row to unit
% numdis = numdis./repmat(sqrt(sum(numdis.^2,2)),1,size(numdis,2));
%  normalize each column to unit -- anzhao lie guiyihua
numdis_se = numdis_se./repmat(sqrt(sum(numdis_se.^2,1)),size(numdis_se,1),1);


numdis_se=(numdis_se/max(max(numdis_se)));
% numdis_se=log(numdis_se);

% cmap = othercolor('RdYlBu4');
% cmap = flip (cmap);
% hold on;colormap(gca,cmap);% colormap

hold on;
conlevel = 100;
contourf(0:1:24,0:2:80,numdis_se,conlevel,'Linestyle','none');
clen=100;cmap = nph_saturate(cbrew('nph_alt_jet',clen),1.0);
colormap(gca,cmap);


% % % % % % % % % % % % % % % % % close all;
% % % % % % % % % % % % % % % % % localtime_se= localtime(id_se);% get Local time (24h)
% % % % % % % % % % % % % % % % % speed_se=speed(id_se)./1000;
% % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % speedpeak = zeros(25,1);
% % % % % % % % % % % % % % % % i=1;
% % % % % % % % % % % % % % % % for lt=0:1:24
% % % % % % % % % % % % % % % %         H = histogram(speed_se(find(localtime_se>lt-0.5 ...
% % % % % % % % % % % % % % % %          & localtime_se<=lt+0.5)),20:0.1:40,'Visible','off');
% % % % % % % % % % % % % % % %     %           it will show a figure window, don't be worried
% % % % % % % % % % % % % % % %         f= fit((20.1:0.1:40).',H.Values.','gauss1', ...% a1, b1, c1
% % % % % % % % % % % % % % % %                 'Lower',[1, 21.0, 10], ...
% % % % % % % % % % % % % % % %                 'Start',[1, 28.0, 12], ... 
% % % % % % % % % % % % % % % %                 'Upper',[1, 32.0, 14]);
% % % % % % % % % % % % % % % %         speedpeak(i)=f.b1;
% % % % % % % % % % % % % % % %         i=i+1;
% % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % hold on;plot(0:1:24,speedpeak,'Color','k','linest','--','LineWidth',.08*fs);

%% only add median of speed
speedmed = zeros(25,1);
% speedmean = zeros(25,1);
i=1;
for lt=0:1:24
        Hinds = speed_se(find(localtime_se>lt-0.5 ...
         & localtime_se<=lt+0.5));
        speedmed(i) = median(Hinds,"omitnan");
%         speedmean(i) = nanmean(Hinds);
        i=i+1;
end


hold on;plot(0:1:24,speedmed,'Color','k','LineWidth',.23*fs);
% hold on;plot(0:1:24,speedmean,'Color','w','LineWidth',.23*fs);

sppeak_se = speedmed;
pkmx_se=num2str(round(max(sppeak_se),4,"significant"));
disp(['max = ' pkmx_se ' km/s']);
nph_text([0.63 0.9],['max = ' pkmx_se ' km/s'], 'FontSize',0.9*fs); % ,'textborder','w');

pkmi_se=num2str(round(min(sppeak_se),4,"significant"));
disp(['min = ' pkmi_se ' km/s']);
nph_text([0.63 0.8],['min = ' pkmi_se ' km/s'], 'FontSize',0.9*fs); % ,'textborder','w');

clim([0 1]);%%%% LIMITS
xlim([0 24])
xtick(0:8:24)
xminortick(0:2:24)

ylim([0 80])
ytick(0:20:80)
yminortick(0:4:80)

%%%% LABELS
% ylabel('Speed (km/s)')
xlabel('Local time (hour)')
box on;

end

%% mch_localtime_contourf
function numdis_se = mch_localtime_contourf(id_se, localtime, mch, fs)
dashcolor = '#a5a5a5';

% mcest_se=mcest(id_se);
localtime_se= localtime(id_se);% get Local time (24h)
mch_se=mch(id_se);

numdis_se = zeros(41,25);
n=1;
for ht=70:110
    m=1;
    for hr=0:1:24
        numdis_se(n,m)=length(find(localtime_se>hr-0.5 ...
            & localtime_se<=hr+0.5 & mch_se>ht-0.5 & mch_se<ht+0.5));
%         mcount_se(m)=length(find(localtime_se>hr-0.5 ...
%             & localtime_se<=hr+0.5));
        m=m+1;
    end
    n=n+1;
end

% %  normalize each row to unit
% numdis = numdis./repmat(sqrt(sum(numdis.^2,2)),1,size(numdis,2));
%  normalize each column to unit -- anzhao lie guiyihua
numdis_se = numdis_se./repmat(sqrt(sum(numdis_se.^2,1)),size(numdis_se,1),1);


numdis_se=(numdis_se/max(max(numdis_se)));
% numdis_se=log(numdis_se);

% cmap = othercolor('PRGn4');
% cmap = flip (cmap);
% hold on;colormap(gca,cmap);% colormap

hold on;
conlevel = 100;
contourf(0:1:24,70:110,numdis_se,conlevel,'Linestyle','none');
clen=100;cmap = nph_saturate(cbrew('nph_alt_jet',clen),1.0);
colormap(gca,cmap);


mchpeak = zeros(25,1);
mchwidth = zeros(25,1);
i=1;
for lt=0:1:24
        H = histogram(mch_se(find(localtime_se>lt-0.5 ...
         & localtime_se<=lt+0.5)),70:0.1:110,'Visible','off');
    %           it will show a figure window, don't be worried
        f= fit((70.1:0.1:110).',H.Values.','gauss1');
        mchpeak(i)=f.b1;
        mchwidth(i)=f.c1;
        i=i+1;
end

hold on;plot(0:1:24,mchpeak+mchwidth,'Color','k','linest','--','LineWidth',.08*fs);
hold on;plot(0:1:24,mchpeak,'Color','k','LineWidth',.2*fs);
hold on;plot(0:1:24,mchpeak-mchwidth,'Color','k','linest','--','LineWidth',.08*fs);

clim([0 1]);%%%% LIMITS
xlim([0 24])
xtick(0:8:24)
xminortick(0:2:24)

ylim([70 110])
ytick(70:10:110)
yminortick(70:2:110)

% % % % add dashed lines to the max hours
% % maxmc=round(maxmc,2,"significant");
% % hold on; plot([maxmc maxmc],[70 170],'linest','--','color','k','LineWidth',.08*fs);
% % maxc=num2str(maxmc);
% % disp(['t_{peak}=' maxc 'hr']);
% % nph_text([0.63 0.80],['peakhour=' maxc 'hr'], 'FontSize',0.9*fs); % ,'textborder','w');

%%%% LABELS
[mcm_se,mcn_se]=hist(mch_se,70:1:110);
[htwidth_se, htpeak_se, ~]=mygaussfit(mcn_se,mcm_se);

% mnum=num2str(round(metnum_se));
% nph_text([0.63 0.87],['n=' mnum], 'FontSize',0.9*fs); % ,'textborder','w');

% hold on; plot([0 24],[htpeak_se htpeak_se],'linest','--','color','k','LineWidth',.08*fs);
% htadd_se=htpeak_se+htwidth_se;
% hold on; plot([0 24],[htadd_se htadd_se],'linest','--','color','k','LineWidth',.08*fs);
% htsub_se=htpeak_se-htwidth_se;
% hold on; plot([0 24],[htsub_se htsub_se],'linest','--','color','k','LineWidth',.08*fs);

hpk_se=num2str(round(htpeak_se,4,"significant"));
disp(['\mu = ' hpk_se ' km']);
nph_text([0.63 0.8],['\mu = ' hpk_se ' km'], 'FontSize',0.9*fs); % ,'textborder','w');

hwt_se=num2str(round(htwidth_se,4,"significant"));
disp(['\sigma = ' hwt_se ' km']);
nph_text([0.63 0.9],['\sigma = ' hwt_se ' km'], 'FontSize',0.9*fs); % ,'textborder','w');

box on;

end

%% mcount_localtime_hist
function  mcount_localtime_hist(id_se, localtime)
dashcolor = '#a5a5a5';
grid on;
% mcest_se=mcest(id_se);
localtime_se= localtime(id_se);% get Local time (24h)
% localtime_se= mod(mcest_se+8,24);% get Local time (24h)

%figure;
hold on;
H = histogram(localtime_se,0:0.1:24,'visible','off');
% H = histogram(localtime_se,0:0.1:24,'visible','on','color','');
binedges = H.BinEdges(1:end-1);
H.BinCounts = H.BinCounts/100;

clen=24;
% plot by time of day:
csteps = 0:(24/clen):(24-(24/clen));
% cmap = nph_saturate(cbrew('nph_Parula',clen),1.0);
cmap = nph_saturate(cbrew('nph_CyclicRainbow',clen),0.8);

% colour by time:
for i = 1:clen
    cinds = inrange(binedges,[csteps(i) csteps(i)+(24/clen)]);
    hold on; b=bar(binedges(cinds),H.Values(cinds), ...
        'facecolor',cmap(clen-i+1,:),'Linestyle','none');
    b.EdgeColor = 'none';
end

%%%% LIMITS
xlim([0 24])
xtick(0:8:24)
xminortick(0:2:24)
xlabel('Local time (hour)')

%ylim([70 110])
%ytick(70:10:110)
%yminortick(70:2:110)

% add dashed lines to the max hour
[ylim1,maxmc] = max(H.BinCounts);
maxmc=0+0.1*maxmc;
hold on; plot([maxmc maxmc],[0 ylim1],'linest','--','color','k','LineWidth',.08*fs);
maxmc=round(maxmc,2,"significant");
maxc=num2str(maxmc);
disp(['t_{peak}=' maxc ' hr']);
nph_text([0.63 0.75],['t_{peak}=' maxc ' hr'], 'FontSize',0.9*fs); % ,'textborder','w');

metnum_se=length(find(id_se));
mnum=num2str(round(metnum_se));
disp(['Meteor Number=' mnum]);
nph_text([0.63 0.88],['n=' mnum], 'FontSize',0.9*fs); % ,'textborder','w');

box on;
end
