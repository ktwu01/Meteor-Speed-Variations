


% CreateOptionalData_V4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ AND TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load original 2014-2022 data
clear all;clc;close all;
load('D:\0lrn\00Res\Data\MMRdata_2014_2022.mat');

addpath 'D:\0lrn\00Res\Functions'

% % % % % % num: start from 2014-01-01
% % % % % % est: start from 2014-01-01 0:00 UT
% % % % % % recorded meteors: start from 74rd day 3 hour of 2014 ???
% % % % % % 24*(datenum(2022,12,31)- datenum(2014,01,01)+1)=max(mcest) 
% % % % % % min(mcest)=1.7575e+03 (about 73 days * 24 hours=1572)

% clear all the abnormal data
tmp1=find(mcmet.mspeed<3000);
        mcmet.est(tmp1)=nan;
        mcmet.hmet(tmp1)=nan;
        mcmet.da(tmp1)=nan;
        mcmet.zenith(tmp1)=nan;
        mcmet.azimuth(tmp1)=nan;
        mcmet.mspeed(tmp1)=nan;
        mcmet.range(tmp1)=nan;
        mcmet.dt(tmp1)=nan;
        mcmet.rv(tmp1)=nan;
        
% define data for all operational days
% BUT use all number to histogram number
%%%%%%% tmp= find(mcmet.est>24*(dayrange(1)- datenum(2014,01,01)) & mcmet.est<24*(dayrange(end)- datenum(2014,01,01)))
        mcest=mcmet.est;
        mch=mcmet.hmet;
%         mcda=mcmet.da;
        mczenith=mcmet.zenith;
        azimuth=mcmet.azimuth;
        speed=mcmet.mspeed;
        mcr=mcmet.range;
%         mcdt=mcmet.dt;
%         mcrv=mcmet.rv;
        number=mcmet.number;

        mcz=pi*mczenith/180;
        mcrd=mcr.*sin(mcz); %groud range
        
        localtime = mod(mcest+8,24); % Beijing time
%         hourofUT = mod(floor(mcest),24);

        dynum = floor(mcest/24)+datenum(2014,01,01); % 0 - 3286 + datenum(2014,01,01)
        doy = daynumber(dynum); % use daynumber function in NPH functions
        % doy: 1-366
        
        % [nmczenith,mmczenith]=hist(mczenith,0:3:90);
        % [nmcazimuth,mmczimuth]=hist(mcazimuth,5:15:360);
        % [mcx,mcy]=pol2cart(mcazimuth,mcrd);
        % rdmc=sqrt(mcx.^2+mcy.^2);
        % [nnmc,mmmc]=hist(rdmc,5:15:360);
        % [nnum,mnum]=hist(mcest,0.5:1:143.5);

        % calculate the days when the MMR is in operation
%         ndays = datenum(2022,12,31)-datenum(2014,01,01)+1-73;

clear mcmet;

%% load original 2023 data

load('D:\0lrn\00Res\Data\Mengcheng_meteor_data_2023.mat');

% % % % % % num: start from 2023-01-01
% % % % % % est: start from 2023-01-01 0:00 UT

% clear all the abnormal data
tmp1=find(mcmet.mspeed<3000);
        mcmet.est(tmp1)=nan;
        mcmet.hmet(tmp1)=nan;
        mcmet.da(tmp1)=nan;
        mcmet.zenith(tmp1)=nan;
        mcmet.azimuth(tmp1)=nan;
        mcmet.mspeed(tmp1)=nan;
        mcmet.range(tmp1)=nan;
        mcmet.dt(tmp1)=nan;
        mcmet.rv(tmp1)=nan;
        
% define data for all operational days
% BUT use all number to histogram number
%%%%%%% tmp= find(mcmet.est>24*(dayrange(1)- datenum(2023,01,01)) & mcmet.est<24*(dayrange(end)- datenum(2023,01,01)))
        mcest1=mcmet.est + 24*(datenum(2023,1,1)- datenum(2014,1,1));

figure; plot(mcest); disp(max(mcest));
figure; plot(mcest1); disp(min(mcest1));

        mch1=mcmet.hmet;
%         mcda1=mcmet.da;
        mczenith1=mcmet.zenith;
        azimuth1=mcmet.azimuth;
	speed1=mcmet.mspeed;
        mcr1=mcmet.range;
%         mcdt1=mcmet.dt;
%         mcrv1=mcmet.rv;
        number1=mcmet.number;

        mcz1=pi*mczenith1/180;
        mcrd1=mcr1.*sin(mcz1); %groud range
        
        localtime1 = mod(mcest1+8,24); % Beijing time
%         hourofUT1 = mod(floor(mcest1),24);

        dynum1 = floor(mcest1/24)+datenum(2023,01,01); % 0 - 3286 + datenum(2023,01,01)
        doy1 = daynumber(dynum1); % use daynumber function in NPH functions

%% Connect
        mcest = [mcest, mcest1];
        mch = [mch, mch1];
        mczenith = [mczenith, mczenith1];
        azimuth = [azimuth, azimuth1];
        speed = [speed, speed1];
        mcr = [mcr, mcr1];
        number = [number, number1];
        mcrd = [mcrd, mcrd1];
        localtime = [localtime, localtime1];
        doy = [doy, doy1];

%% Transpose, so that they are column vectors
        mcest = mcest';
        mch = mch';
        mczenith = mczenith';
        azimuth = azimuth';
        speed = speed';
        mcr = mcr';
        number = number';
        mcrd = mcrd';
        localtime = localtime';
        doy = doy';

save('D:\0lrn\00Res\Data\McMRdata_2014_2023.mat', ...
    'mcest','mch','mczenith', ...
    'azimuth','speed','mcr','number','mcrd','localtime','doy');

disp('All Data Saved! Please Check!')

% % filename=[''];
% % MH_Zonal=[];MH_Meridional=[]; %% 漠河定义风场数组 dingyi array

% % % % Apply Zenith Angle Limits:
% % % MWD.Thresholds.ZenithLimits = [15 65]; % use [15 65] for computing mean winds, impose stricter for GWs.
% % % inds = zen >= MWD.Thresholds.ZenithLimits(1) & zen <= MWD.Thresholds.ZenithLimits(2);
