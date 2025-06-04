%% Earth rotation speed around the NP Use SI
clear all;close all
% calculation

lat = 33.4 * pi/ 180; % rad
R = 6371e3; %Earth radius, m
alt = 0; %altitude, m
T = 24*3600; %Earth rotation period, s
spe_rot = cos(lat)*(R+alt)*2*pi/T; % m/s
spe_rotkms = spe_rot/1000;

disp(['Earth rotation speed at ',num2str(lat), ' latitude is ',num2str(spe_rot) ' m/s']);
disp(['Earth rotation speed at ',num2str(lat), ' latitude is ',num2str(spe_rotkms) ' km/s'])