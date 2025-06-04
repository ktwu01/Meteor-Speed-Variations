




%%% WKT_RgHtSp_4col_V1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ AND TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % num: start from 2014-01-01
% % % % % % est: start from 2014-01-01 0:00 UT

% % % % % % 24*(datenum(2022,12,31)- datenum(2014,01,01)+1)=max(mcest) 
% % % % % % min(mcest)=1.7575e+03 (73 days * 24 hours=1752)

% % % % % % recorded meteors: start from 74rd day 5.5 hour of 2014

clear all;clc;close all; time_sys3=clock;
% load('../../Data/OptionalMMRdata2014_2022.mat');
% load(['../../Data/ids2014_2022.mat']);
load('../../Data/MCMR_MetData_2014_2023.mat');
load('../../Data/MCMR_ids_2022.mat');

addpath '../../Functions'

figure; hold all;
whitefig;

%--------------------------------------------------------
rows = 3; cols = 4;

figpos([.9 1]);
vert_gap = 0.045*sqrt(5/rows);    horz_gap = 0.02*sqrt(10/cols);
lower_marg = 0.05*sqrt(5/rows);     upper_marg = 0.02*sqrt(5/rows)*1.8;
left_marg = 0.035*sqrt(10/cols);     right_marg = 0.04*sqrt(10/cols);
subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);
%--------------------------------------------------------
fs = 20; % fontsize

%% barcolor bar color
% Example: barcolor = '#0065a1';
dashcolor = '#a5a5a5';
barcolor = '#009BCA';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SubPlots : 1-2 rols contourf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% id_pr : inds for Perihelion Jan-4-2022 R_earthandsun= 147105052 km
% id_ap : inds for Aphelion July-4-2022 R_earthandsun= 152098455 km

%% plot
subplot(rows,cols, 1);numdis_se1 = mcr_speed_contourf(id_pr, speed, mcr);ylabel('Range (km)');title('Jan 4');% clim([0 1]);
subplot(rows,cols, 2);numdis_se2 = mcr_speed_contourf(id_ap, speed, mcr);title('Apr 4');% clim([0 1]);
subplot(rows,cols, 3);numdis_se3 = mcr_speed_contourf(id_pr, speed, mcr);title('Jul 4');% clim([0 1]);
subplot(rows,cols, 4);numdis_se4 = mcr_speed_contourf(id_ap, speed, mcr);title('Oct 2');% clim([0 1]);

% %%%% COLORBAR
% hold on;
% cbar = nph_colorbar;clim([0 1]);
% cbar.Label.String = ['N/N_{max}']; %cbar.Label.String = ['log_{10}(N/N_{max})'];
% cbar.TickDirection = 'out';
% cbar.Position = cbar.Position .* [1 1 0.8 1];
% cbar.Position = cbar.Position + [0.5*cbar.Position(3) 0 0 0];
% cbar.Label.Rotation = 90;

%% plot
subplot(rows,cols, 5);numdis_se5 = mch_speed_contourf(id_pr, speed, mch, fs);ylabel('Altitude (km)');% clim([0 2000]);
subplot(rows,cols, 6);numdis_se6 = mch_speed_contourf(id_se, speed, mch, fs);% clim([0 2000]);
subplot(rows,cols, 7);numdis_se7 = mch_speed_contourf(id_ap, speed, mch, fs);% clim([0 2000]);
subplot(rows,cols, 8);numdis_se8 = mch_speed_contourf(id_ae, speed, mch, fs);% clim([0 2000]);

%%%% COLORBAR
hold on;
cbar = nph_colorbar;clim([0 1]);
cbar.Label.String = ['N/N_{max}']; %cbar.Label.String = ['log_{10}(N/N_{max})'];
cbar.TickDirection = 'out';
cbar.Position = cbar.Position .* [1 1 0.8 1];
cbar.Position = cbar.Position + [0.5*cbar.Position(3) 14*cbar.Position(3) 0 0];
%             cbar.Label.String = ['\bf{' clabs{ax} '}'];
cbar.Label.Rotation = 90;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Meteors speed distribution (9 yrs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

len = 80;
maxsp = 80;
minsp = 0;

% ax = 9;
subplot(rows,cols, 9);y1 = sphist(speed,id_pr,len,maxsp,minsp, fs);%title('Jan 4');
ylabel('N/N_{max}')
% hold on; nph_text([0.015 0.875],['(' alphabet(ax) ')'],'fontsize',1.5*fs); 
% ax = ax + 1; 

subplot(rows,cols, 10);y2 = sphist(speed,id_se,len,maxsp,minsp, fs);%title('Apr 4');
% hold on; nph_text([0.015 0.875],['(' alphabet(ax) ')'],'fontsize',1.5*fs); 
% ax = ax + 1; 

subplot(rows,cols, 11);y3 = sphist(speed,id_ap,len,maxsp,minsp, fs);%title('Jul 4');
% hold on; nph_text([0.015 0.875],['(' alphabet(ax) ')'],'fontsize',1.5*fs); 
% ax = ax + 1; 

subplot(rows,cols, 12);y4 = sphist(speed,id_ae,len,maxsp,minsp, fs);%title('Oct 2');
% hold on; nph_text([0.015 0.875],['(' alphabet(ax) ')'],'fontsize',1.5*fs); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GLOBAL FORMATTING FOR ALL AXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ax = 1:rows*cols
%     subplot(rows,cols,axs{ax});
    subplot(rows,cols,ax);

    axx = gca;
    setfont(fs);
    set(gca,'linewi',1,'tickdir','out');

hold on; nph_text([0.015 0.875],['(' alphabet(ax) ')'],'fontsize',1.5*fs); 
end

time_sys4=clock;

timerun_43=etime(time_sys4,time_sys3);
time43=num2str(floor(timerun_43));

error

%% EXPORT FIG ==============================================================
% Set renderer to "painters"
set(gcf, 'Renderer', 'painters');
figure_name = ['Fig8'];
print(gcf,figure_name,'-djpeg','-r600');
% print(gcf,['WKT_RgHtSp_4col_V1'],'-dpng','-r600');
% saveas(gcf,'Figure8','epsc')
% saveas(gcf,['WKT_RgHtSp_4col_V1_t=',time43,'s'],'svg');
disp('Figure Saved.');

%% mcr_speed_contourf
function numdis_se = mcr_speed_contourf(id_se, speed, mcr)

mcr_se=mcr(id_se);
speed_se=speed(id_se)./1000;

numdis_se = zeros(42,41);
n=1;
for rg=70:6:320
    m=1;
    for sp=0:2:80
        numdis_se(n,m)=length(find(speed_se>sp-1 ...
         & speed_se<=sp+1 & mcr_se>rg-3 & mcr_se<rg+3));
        m=m+1;
    end
    n=n+1;
end

% numdis_se=numdis_se'
% % numdis_se=(numdis_se/max(max(numdis_se)));

hold on;
conlevel = 100;
contourf(0:2:80,70:6:320,numdis_se,conlevel,'Linestyle','none');
clen=100;cmap = nph_saturate(cbrew('nph_alt_jet',clen),1.0);
colormap(gca,cmap);

%%%% LIMITS
xlim([0 80])
xtick(0:20:80)
xminortick(0:5:80)

ylim([70 310])
ytick(70:60:310)
yminortick(70:10:310)

box on; %%%% LABELS
% ylabel('Range (km)')

end

%% mch_speed_contourf
function numdis_se = mch_speed_contourf(id_se, speed, mch, fs)

mch_se=mch(id_se);
speed_se=speed(id_se)./1000;

numdis_se = zeros(41,41);
n=1;
for ht=70:110
    m=1;
    for sp=0:2:80
        numdis_se(n,m)=length(find(speed_se>sp-1 ...
         & speed_se<=sp+1 & mch_se>ht-0.5 & mch_se<ht+0.5));
        m=m+1;
    end
    n=n+1;
end

% numdis_se=numdis_se'
numdis_se=(numdis_se/max(max(numdis_se)));

hold on;conlevel = 100;
contourf(0:2:80,70:110,numdis_se,conlevel,'Linestyle','none');
hold on;nph_text([0.72 0.05],['# = ' num2str(length(id_se))],'fontsize',fs*0.9);
clen=100;cmap = nph_saturate(cbrew('nph_alt_jet',clen),1.0);
colormap(gca,cmap);
clim([0 1]);

% % % % mchpeak = zeros(41,1);
% % % % i=1;
% % % % for sp=0:2:80
% % % % %     [mcm,mcn]=hist(mch_se(find(speed_se>sp-1 ...
% % % % %             & speed_se<=sp+1)),70:0.1:110);
% % % % %     [~, mchpeak(i), ~]=mygaussfit(mcn,mcm);
% % % %         H = histogram(mch_se(find(speed_se>sp-1 ...
% % % %          & speed_se<=sp+1)),70:0.1:110,'Visible','off');
% % % %     %           it will show a figure window, don't be worried
% % % %         f= fit((70.1:0.1:110).',H.Values.','gauss1');
% % % %         mchpeak(i)=f.b1;
% % % %         i=i+1;
% % % % end
% % % % 
% % % % % find the real part of height peak and width
% % % % mchpeak=real(mchpeak);
% % % % % mchwidth=real(mchwidth);
% % % % hold on;plot(12:2:70,mchpeak(6:1:35),'Color','k','LineWidth',2);
% % % % 
% % % % % mednum = nanmedian(numdis_se,1);
% % % % % hold on;plot(0:2:80,mednum,'Color','k','LineWidth',2);

%%%% LIMITS
xlim([0 80])
xtick(0:20:80)
xminortick(0:5:80)

ylim([70 110])
ytick(70:10:110)
yminortick(70:2:110)

box on;
%%%% LABELS
% ylabel('Altitude (km)')

end



function y = sphist(speed,id_se,len,maxsp,minsp, fs)
    barcolor = '#fbe3a4';

    speed_se = speed(id_se);
    H = histogram(speed_se./1000,minsp:(maxsp-minsp)/len:maxsp,'visible','off');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% fit
    % step1: normalize the number of meteors observed in a day
    % so that we can control the parameters in an acceptable region
    y = H.BinCounts/max(H.BinCounts);
    x = H.BinEdges(1:end-1);

    bar(x,y,'barwidth',1,'facecolor',barcolor)
    %%%% LIMITS
    xlim([minsp maxsp])
    xtick(minsp:20:maxsp)
    xminortick(minsp:5:maxsp)
    xlabel('Speed (km/s)');

    ylim([0 1.1])
    ytick(0:0.2:1)
    yminortick(0:0.1:1)

            %%%%%%%%%%%%%%%%%%%%%%%%%% fit
    f = fit(x.',y.','gauss2', ...% a1, mu1, sig1, a2, mu2, sig2
                    'Lower',[1, 24.0, 5, 0.10, 44, 8], ...
                    'Start',[1, 28.0, 10, 0.40, 54, 11], ...
                    'Upper',[1, 40.0, 15, 0.46, 60, 14]);

    yi= f.a1*exp(-((x-f.b1)/f.c1).^2) + f.a2*exp(-((x-f.b2)/f.c2).^2);
    hold on;plot(x,yi,'linest','--','color','r','LineWidth',.2*fs);

    disp(['Gauss2 fit result']);
    disp(['a1=' num2str(round(f.a1,4,"significant"))]);
%     nph_text([0.72 0.85],['a_1=' num2str(round(f.a1,4,"significant"))], 'FontSize',0.9*fs); %,'textborder','w');
    disp(['a2=' num2str(round(f.a2,4,"significant"))]);
    nph_text([0.72 0.84],['a_2 = ' num2str(round(f.a2,4,"significant"))], 'FontSize',0.9*fs); %,'textborder','w');

    sppk1=num2str(round(f.b1,4,"significant"));
    disp(['b1=' sppk1 ]);
    nph_text([0.74 0.74],['\mu_{1} = ' sppk1 ], 'FontSize',0.9*fs); %,'textborder','w');
    
    sppk2=num2str(round(f.b2,4,"significant"));
    disp(['b2=' sppk2 ]);
    nph_text([0.74 0.64],['\mu_{2} = ' sppk2 ], 'FontSize',0.9*fs); %,'textborder','w');

    spwt1=num2str(round(f.c1,4,"significant"));
    disp(['c1=' spwt1 ]);
    % disp(['\sigma_{1}=' spwt1 ]);
    nph_text([0.74 0.54],['\sigma_1 = ' spwt1  ], 'FontSize',0.9*fs); %,'textborder','w');

    spwt2=num2str(round(f.c2,4,"significant"));
    disp(['c2=' spwt2 ]);
    nph_text([0.74 0.44],['\sigma_2 = ' spwt2 ], 'FontSize',0.9*fs); %,'textborder','w');

end