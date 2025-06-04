



% Calcu_Speed_Gauss2_ave_V4
% load original data
clear all;clc;close all; time_sys1=clock;
load('D:\000learning\00research\00DataAndCode\Mengcheng_wkt_Nhindley\OptionalMMRdata2014_2022.mat');

% Mention! BEFORE YOU RUN THIS CODE
% please turn on ALL or SOME of the if (1) so that you can satisfile your needs!
% try not to just directely run this file instead.

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
            H=histogram(speed(inds)./1000,0:0.1:clen/10,'visible','off');
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
save('D:\000learning\00research\00DataAndCode\Mengcheng_wkt_Nhindley\mcspgs2_ave_2014_2022.mat', ...
    'spa1ave','spa2ave','sppk1ave','sppk2ave',"spwt1ave","spwt2ave");
disp('speed peak and width Data Saved! Please Check!')

end

%% calcu the time of this code
time_sys2=clock;

timerun_21=etime(time_sys2,time_sys1);
time21=num2str(floor(timerun_21));
disp(['time21=',time21,'s'])