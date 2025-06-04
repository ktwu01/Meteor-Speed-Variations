% load original data
clear all;clc;close all;
load('D:\0lrn\00Res\Data\McMR_WKT_Data_2014_2023.mat','mch','mcest');

% Mention! BEFORE YOU RUN THIS CODE
% please turn on ALL or SOME of the if (1) so that you can satisfile your needs!
% try not to just directely run this file instead.
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %% height contourf of each month of years
% % % % time use: about 7 mins
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % if(0)
% % %     n=1;
% % % for ht=70:110
% % %     m=1;
% % %     for mt=1:108 %month 12*9=108
% % %         htdis9yrs(n,m)=length(find(mcest>720*(mt-1) & mcest<=720*(mt) ...
% % %             & mch>ht-0.5 & mch<ht+0.5));
% % %         m=m+1;
% % %     end
% % %     n=n+1;
% % % end
% % % 
% % % disp('Caculated: height contourf data of each month of years')
% % % 
% % % % save the file
% % % save('D:\000learning\00research\00DataAndCode\Mengcheng_wkt_Nhindley\mchdis_2014_2022.mat', ...
% % %     'htdis9yrs');
% % % 
% % % disp('Saved: height contourf data of each month of years.')
% % % disp('Please Check whether it has been truly saved!')
% % % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate the height peak and Width of each day of years 
% time use: about 6-7 mins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1)
    DAY_PER = datenum(2023,12,31)-datenum(2014,1,1)+1;

mchpeak=zeros(DAY_PER,1);
mchwidth=zeros(DAY_PER,1);
    for i = 1:DAY_PER
        z = find(mcest >=24*(i-1) & mcest <24*i);
        if z>1500 % more than 1500 per day will be included
% % % % %             mcest_daily=mcest (z);
            h_daily=mch(z);
%             [mcm,mcn]=hist(mch(z),70:0.4:110);
%             [mchwidth(i), mchpeak(i), htA]=mygaussfit(mcn,mcm);
% % % % %           histogram
% % % % %           H = histogram(h_daily,60:1:120);
            H = histogram(h_daily,70:0.4:110,'Visible','off');
% % % % %           it will show a figure window, don't be worried
% % % % %           fit
            f= fit((70.4:0.4:110).',H.Values.','gauss1');
            mchpeak(i)=f.b1;
            mchwidth(i)=f.c1;
% % % %             mchwidth(i)=f.c1/sqrt(2);
        else % less than 1500 per day will be excluded
            mchwidth(i)=nan;
            mchpeak(i)=nan;
        end
    end

% find the real part of height peak and width
mchpeak=real(mchpeak);
mchwidth=real(mchwidth)/sqrt(2);
disp('calculate the height peak and Width of each day of years Done!')

% save the file
save('D:\0lrn\00Res\Data\Mch_pk_wid_gs1_2014_2023.mat', ...
    'mchwidth','mchpeak');
disp('height peak and width Data Saved! Please Check!')

end