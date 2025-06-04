



% Calcu_Speed_Gauss2_3287.m
% load original data
clear all;clc;close all; time_sys1=clock;
load('D:\000learning\00research\00DataAndCode\Mengcheng_wkt_Nhindley\OptionalMMRdata2014_2022.mat');

% Mention! BEFORE YOU RUN THIS CODE
% please turn on ALL or SOME of the if (1) so that you can satisfile your needs!
% try not to just directely run this file instead.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate the speed peak and Width of each day of 9 years 
% time use: about 500s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)
            spa2=zeros(3287,1);
            spa1=zeros(3287,1);            
            spwt2=zeros(3287,1);
            sppk2=zeros(3287,1);            
            spwt1=zeros(3287,1);
            sppk1=zeros(3287,1);

    for i = 1:3287 % Caution! there isn't any data during Jan- Mar 2014
        inds = find(mcest >=24*(i-1) & mcest <24*i);
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

figure;scatter(1:3287,spa2);xlim([0 3287]);title('spa2');
print('spa2.png', '-dpng')

figure;scatter(1:3287,spa1);xlim([0 3287]);title('spa1');
print('spa1.png', '-dpng')

figure;scatter(1:3287,spwt2);xlim([0 3287]);title('spwt2');
print('spwt2.png', '-dpng')

figure;scatter(1:3287,sppk2);xlim([0 3287]);title('sppk2');
print('sppk2.png', '-dpng')

figure;scatter(1:3287,spwt1);xlim([0 3287]);title('spwt1');
print('spwt1.png', '-dpng')

figure;scatter(1:3287,sppk1);xlim([0 3287]);title('sppk1');
print('sppk1.png', '-dpng')

disp('calculate the speed peak and Width of each day of 9 years Done!')

% save the file
save('D:\000learning\00research\00DataAndCode\Mengcheng_wkt_Nhindley\mcspgs2_9yrs_2014_2022.mat', ...
    'spa1','spa2','sppk1','sppk2',"spwt1","spwt2");
disp('speed peak and width Data Saved! Please Check!')

end

%% calcu the time of this code

time_sys2=clock;

timerun_21=etime(time_sys2,time_sys1);
time21=num2str(floor(timerun_21));
disp(['time21=',time21,'s'])
