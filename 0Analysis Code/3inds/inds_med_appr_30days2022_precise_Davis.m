



%% calculate specific days

% ecc
e = 0.0167;
a = 1;
c = a*e;
b = sqrt(a*a-c*c);

% theta: 0 , nishizhen from perihelion (Jan 4)
% d t(theta) = d s(theta) / v(theta)
dt = @(theta) cos(theta).*cos(theta)+(1-e*e)*sin(theta).*sin(theta);


% calcu: from perihelion to spring
theta1  = atan (b/c)
delta_t_1 = integral(dt,0,theta1)
% calcu: from spring to aphelion
delta_t_2 = integral(dt,theta1,pi)

% calcu: proportion to  1 year time
delta_time1 = 365.256*0.5*delta_t_1/(delta_t_1+delta_t_2)
delta_time2 = 365.256*0.5*delta_t_1/(delta_t_1+delta_t_2)


% calcu: proportion to  1 year time
perihl_datenum = datenum(2005,01,04)-1
spring_datenum = datestr((datenum(2005,01,04)+delta_time1))
autumn_datenum = datestr((datenum(2005,07,04)+delta_time2))

%% calculate: find +- 15 days around those specific days
load('D:\0lrn\00Res\Data\davismet_dataonly2005.mat','d33met');
% % % id_pr : inds for Perihelion Jan-4-2022 R_earthandsun= 147105052 km
% % % id_ap : inds for Aphelion July-4-2022 R_earthandsun= 152098455 km
mcest = d33met.est;
if(1)
id_pr = find(mcest>24*(datenum(2005,01,04)- datenum(2005,01,01)-15) ...
    & mcest<24*(datenum(2005,01,04)- datenum(2005,01,01)+15));

id_se = find(mcest>24*(datenum(spring_datenum)- datenum(2005,01,01)-15) ...
    & mcest<24*(datenum(spring_datenum)- datenum(2005,01,01)+15));
    
id_ap = find(mcest>24*(datenum(2005,07,04)- datenum(2005,01,01)-15) ...
    & mcest<24*(datenum(2005,07,04)- datenum(2005,01,01)+15));

id_ae = find(mcest>24*(datenum(autumn_datenum)- datenum(2005,01,01)-15) ...
    & mcest<24*(datenum(autumn_datenum)- datenum(2005,01,01)+15));

save(['D:\0lrn\00Res\Data\ids_Davis_2005.mat'], ...
    'id_se','id_ae','id_pr','id_ap' ...
    );
disp ('saved');
end