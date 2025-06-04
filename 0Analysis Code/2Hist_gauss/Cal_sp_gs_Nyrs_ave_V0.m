



% Calcu_Speed_Gauss2_DAY_PER.m
% load original data
clear all;clc;close all;
addpath 'D:\0lrn\00Res\Functions'
load('D:\0lrn\00Res\Data\McMRdata_2014_2023.mat','speed','mcest','doy');
speed = speed/1000;
DAY_PER = datenum(2023,12,31)- datenum(2014,1,1)+1;

% Mention! BEFORE YOU RUN THIS CODE
% please turn on ALL or SOME of the if (1) so that you can satisfile your needs!
% try not to just directely run this file instead.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate the speed peak and Width of each day of 9 years 
% time use: about 500s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)
            spa2=zeros(DAY_PER,1);
            spa1=zeros(DAY_PER,1);            
            spwt2=zeros(DAY_PER,1);
            sppk2=zeros(DAY_PER,1);            
            spwt1=zeros(DAY_PER,1);
            sppk1=zeros(DAY_PER,1);

    for i = 1:DAY_PER % Caution! there isn't any data during Jan- Mar 2014
        inds = find(mcest >=24*(i-1) & mcest <24*i);
        if length(inds)>1500 % more than 1500 per day will be included
            clen=800;
            H=histogram(speed(inds),0:1/10:clen/10,'visible','off');
            binedges = H.BinEdges(1:end-1);
            %%%%%%%%%%%%%%%%%%%%%%%%%% fit
            % step1: normalize the number of meteors observed in a day
            % so that we can control the parameters in an acceptable region
            y = H.BinCounts(1:800)/max(H.BinCounts(1:800));
            x = 0:0.1:80-0.1;
            % y = H.BinCounts(100:800);
            % x = 10:0.1:80;
            % A1, b1, c1, A2, b2, c2
            %%%%%%%%%%%%%%%%%%%%%%%%%% fit
    f = fit(x.',y.','gauss2', ...% a1, mu1, sig1, a2, mu2, sig2
                    'Lower',[1, 24.0, 5, 0.10, 44, 8], ...
                    'Start',[1, 28.0, 10, 0.40, 54, 11], ... 
                    'Upper',[1, 40.0, 15, 0.46, 60, 14]);

%     yi= f.a1*exp(-((x-f.b1)/f.c1).^2) + f.a2*exp(-((x-f.b2)/f.c2).^2);

            % % % % % % % % % % % % % % % % % test the fit result
%             figure;bar (x,y);
%             hold on;plot(x,yi,'linest','--','color','r','linewi',2);
            spa2(i)=f.a2;
            spa1(i)=f.a1;            
            spwt2(i)=f.c2;
%             disp(['\sigma_{2}=' spwt2 ' km/s']);
            % nph_text([0.6 0.9],['\sigma_{2}=' spwt2 ' km/s'],'FontSize',0.8*fs,'textborder','w');
            
            sppk2(i)=f.b2;
%             disp(['\mu_{2}=' sppk2 ' km/s']);
            % nph_text([0.6 0.7],['\mu_{2}=' sppk2 ' km/s'],'FontSize',0.8*fs,'textborder','w');
            
            spwt1(i)=f.c1;
%             disp(['\sigma_{1}=' spwt1 ' km/s']);
            % nph_text([0.6 0.11],['\sigma_{1}=' spwt1 ' km/s'],'FontSize',0.8*fs,'textborder','w');
            
            sppk1(i)=f.b1;
%             disp(['\mu_{1}=' sppk1 ' km/s']);
            % nph_text([0.6 0.01],['\mu_{1}=' sppk1 ' km/s'],'FontSize',0.8*fs,'textborder','w');

        else % less than 1500 per day will be excluded
            spa2(i)= nan;
            spa1(i)= nan;            
            spwt2(i)=nan;
            sppk2(i)=nan;            
            spwt1(i)=nan;
            sppk1(i)=nan;
        end
    end

figure;scatter(1:DAY_PER,spa2);xlim([0 DAY_PER]);title('spa2');
print('spa2.png', '-dpng')

figure;scatter(1:DAY_PER,spa1);xlim([0 DAY_PER]);title('spa1');
print('spa1.png', '-dpng')

figure;scatter(1:DAY_PER,spwt2);xlim([0 DAY_PER]);title('spwt2');
print('spwt2.png', '-dpng')

figure;scatter(1:DAY_PER,sppk2);xlim([0 DAY_PER]);title('sppk2');
print('sppk2.png', '-dpng')

figure;scatter(1:DAY_PER,spwt1);xlim([0 DAY_PER]);title('spwt1');
print('spwt1.png', '-dpng')

figure;scatter(1:DAY_PER,sppk1);xlim([0 DAY_PER]);title('sppk1');
print('sppk1.png', '-dpng')

disp('calculate the speed peak and Width of each day of 9 years Done!')

% save the file
save('D:\0lrn\00Res\Data\Mc_sp_GS2_daily_2014_2023.mat', ...
    'spa1','spa2','sppk1','sppk2','spwt1','spwt2');
disp('speed peak and width Data Saved! Please Check!')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate the speed peak and Width of ave year
% time use: about 2 mins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)
            spa2ave=zeros(366,1);
            spa1ave=zeros(366,1); 
            spwt2ave=zeros(366,1);
            sppk2ave=zeros(366,1); 
            spwt1ave=zeros(366,1);
            sppk1ave=zeros(366,1);

    for i = 1:366 % Caution! there isn't any data during Jan- Mar 2014
        inds = find(doy >=(i-0.5) & doy <i+0.5);
        if length(inds)>1500 % more than 1500 per day will be included
            clen=800;
            H=histogram(speed(inds),0:0.1:clen/10,'visible','off');
            binedges = H.BinEdges(1:end-1);
            %%%%%%%%%%%%%%%%%%%%%%%%%% fit
            % step1: normalize the number of meteors observed in a day
            % so that we can control the parameters in an acceptable region
            y = H.BinCounts(1:800)/max(H.BinCounts(1:800));
            x = 0:0.1:80-0.1;

            % y = H.BinCounts(100:800);
            % x = 10:0.1:80;
            % A1, b1, c1, A2, b2, c2
            %%%%%%%%%%%%%%%%%%%%%%%%%% fit
    f = fit(x.',y.','gauss2', ...% a1, mu1, sig1, a2, mu2, sig2
                    'Lower',[1, 24.0, 5, 0.10, 44, 8], ...
                    'Start',[1, 28.0, 10, 0.40, 54, 11], ...
                    'Upper',[1, 40.0, 15, 0.46, 60, 14]);

            spa2ave(i)=f.a2;
            spa1ave(i)=f.a1; 
            spwt2ave(i)=f.c2; 
            sppk2ave(i)=f.b2;
            spwt1ave(i)=f.c1;
            sppk1ave(i)=f.b1;

        else % less than 1500 per day will be excluded
            spa2ave(i)= nan;
            spa1ave(i)= nan; 
            spwt2ave(i)=nan;
            sppk2ave(i)=nan; 
            spwt1ave(i)=nan;
            sppk1ave(i)=nan;
        end
    end

figure;scatter(1:366,spa2ave);xlim([0 366]);title('spa2ave');
figure;scatter(1:366,spa1ave);xlim([0 366]);title('spa1ave');
figure;scatter(1:366,spwt2ave);xlim([0 366]);title('spwt2ave');
figure;scatter(1:366,sppk2ave);xlim([0 366]);title('sppk2ave');
figure;scatter(1:366,spwt1ave);xlim([0 366]);title('spwt1ave');
figure;scatter(1:366,sppk1ave);xlim([0 366]);title('sppk1ave');

% figure;subplot(2,3,1);scatter(1:366,spa1ave);
% hold on; subplot(2,3,2);scatter(1:366,spa2ave);
% hold on; subplot(2,3,3);scatter(1:366,spwt1ave);
% hold on; subplot(2,3,4);scatter(1:366,sppk1ave);
% hold on; subplot(2,3,5);scatter(1:366,spwt2ave);
% hold on; subplot(2,3,6);scatter(1:366,sppk2ave);

disp('calculate the speed peak and Width of ave year Done!')

% save the file
save('D:\0lrn\00Res\Data\Mc_sp_GS2_ave_2014_2023.mat', ...
    'spa1ave','spa2ave','sppk1ave','sppk2ave','spwt1ave','spwt2ave');
disp('speed peak and width Data Saved! Please Check!')

end
