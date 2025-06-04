%% Earth revolution speed around the Sun Use SI
% https://www.zdic.net/ts/fulu/2019/05/33.html
% Basic Astronomy Cheat Sheet: https://zhuanlan.zhihu.com/p/62221622
% In addition to the rotational speed of the Earth spinning on its axis,
% the planet is also speeding at about 29.78 km/s in its revolution
% around the sun once every 365.2425 days.

% The formula for the velocity of an object at some distance r 
% from the Sun is:
% v  =  sqrt[GM*(2/r - 1/a)]
% Where G is the universal gravitational constant, M is the mass of the 
% Sun, and a is the planet's semimajor axis.
% At perihelion, Earth's distance from the Sun is r = a(1-e) 
% and at aphelion, it's r = a(1+e).
% G = 6.673*10-11 N m2/kg2
% M = 1.989*1030 kg
% a = 1.496*1011 m
% e = 0.01671
% So plugging in the numbers, 
% the speed at perihelion is 30,286 km/s and at aphelion it's 29,291 m/s.

clear all;close all
%% calculate the Earth revolution speed
% Jan 4 2022 is perihelion
doy = (1:365)-4;
theta = 2*pi*(doy)/365; % rad. calculate from perihelion

G = 6.67408e-11; %N m2/kg2
M = 1.98855e30; % kg
e = 0.01671;
a = 1.49597871e11; % m
b = a*sqrt(1-e*e);
c = a*e;
p = b*b/c; %jiaozhunju

% wrong!!! we shall calcu the distance from jiaodian
% r = sqrt( (a*cos(theta)).^2 + (b*sin(theta)).^2 );

% radius from focal point
rfp = e*p./(1+e*cos(theta));
spe_revo = sqrt(G*M*(2./rfp - 1/a));

figure;plot(spe_revo)
find(spe_revo > max(spe_revo)-0.0000000001) % max speed in Jan 4
