



%%% zeni_lt_surf_con_1_4p_V0_01

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

subplot(rows,cols,axs{2});az_lt_contourf(id_pr,localtime, azimuth, 1);title('Jan 4')

function az_lt_contourf(id_ae,localtime, azimuth, px)
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
    colormap(gca,cmap);

    cbar= colorbar;clim([0 1]); % nph_colorbar
%     cbar.Label.String = ['N/N_{max}'];
    ylim([0 360]);
    ytick(0:90:360)
    yminortick(0:30:360)

    xlim([0 24]);
    xtick(0:6:24)
    xminortick(0:2:24)

    
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

    hold on;scatter(0:1:24,azi_mn_deg,'Color','k','LineWidth',2);
    % hold on;plot(0:1:24,speedmean,'Color','w','LineWidth',2);

    %%%% LABELS
    xlabel('Local time (hour)');
    ylabel('Azimuth (deg)');
%     mu = azidtdis;

end
