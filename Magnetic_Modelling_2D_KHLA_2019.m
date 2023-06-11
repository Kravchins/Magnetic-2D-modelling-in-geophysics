%% This program is used to create the polygon and output the total magnetic field
%% Any use of this code MUST refer to the original publication:
% Kravchinsky, V.A., Hnatyshin, D., Lysak, B., Alemie, W., 2019. 
% Computation of magnetic anomalies caused by two‐dimensional structures of arbitrary shape: Derivation and Matlab implementation. 
% Geophysical Research Letters, 46(13), 7345-7351, https://doi.org/10.1029/2019GL082767

close all; clear all; clc
% format long

%% Polygon definition
% Polygon = [50,5.86;64.14,20;50,34.14;35.86,20]; % Diamond
% Polygon = [40,10;60,10;60,30;40,30]; % Square
% Polygon = [45,10;55,10;55,20;45,20]; % Another square
% Polygon = [45,10;55,10;60,30;40,30]; % Trapazoid
Polygon = [45,10;55,10;65,20;60,30;50,35;40,30;35,20]; % irregular shape

% Make sure that the polygon vertexes are in the clockwise direction for the
% reference system like it is shown here:
%    
%     -------------> +X
%    |
%    |      B
%    |   A     C
%    |      D
%    |
%    |
%    \/
%   +Z
%  
% where X is directed along a profile and Z is directed to the center of the Earth

%% Distance along profile
    x = 0:1:100; % distance in meters (or km)
    
%% Define inclination, declination, angle from North in degrees

    brsus   = 0.0001;          % Background magnetic susceptibility (MS), SI units
    % For example, sandstone usually has 0.00001<MS<0.025 and basalt has 0.08<MS<0.35 SI units
    % (from Fundamentals of Geophysics, 2020, Lowrie and Fichtner,
    % Cambridge Press, https://www.cambridge.org/core_title/gb/535177)
    bsus    = 0.0008;         % Body MS (SI units)
    inm     = 56000e-9;    % Strength of Earth's field for the area (nT)
    % For example, inm = 56157.4 (nT) = 56157.4e-9 (T); T = kg⋅s−2⋅A−1
 
    mu_0 = 4*pi*1e-7; % Permiability of the free space = 1.25663706 × 10-6 (m⋅kg⋅s−2⋅A−2)
    brm     = brsus*inm/mu_0;      % Background induced magnetization (m⋅A−1)
    % [(unitless⋅T/(m⋅kg⋅s−2⋅A−2) = (kg⋅s−2⋅A−1)/(m⋅kg⋅s−2⋅A−2) = (A⋅m−1)]
%     For a unit volume
%     vol = 1; unit volume in m^3
%     m=vol*bsus*irm/mu_0;  % Magnetic moment of sphere (ignores remnant magnetization effect)
%     magnetization= m/vol;
    Ii      = 70;   % Inclination angle in degrees for induced changes from -90 to +90 degrees
    Di      = 355;  % Declination angle in degrees  for induced, changes from 0 to 360 deg
    C       = 0;    % The angle between + x axis and geographic north. C=0 means the profile directed to the geograph north
    bm      = bsus*inm/mu_0;   % Body induced magnetization (A⋅m−1)
    
    magrem      = 0.01;    % Remnant magnetization (A⋅m−1). 
    % Take a typical value for the rock or better measure a few representative samples in the lab
    % Ir and Dr can be also calculated from a paleomagnetic pole for the study region
    Ir   = -50;     % Remnant Inclination. You may calculate it from the published APWP for different geological periods
    Dr   = 170;       % Remnant Declination. You may calculate it from the published APWP for different geological periods

    dmag    = (bm-brm);   % Magnetic difference between body and background induced magnetizations (A/m)
    Ji = dmag*mu_0*1e9;   % Magnetic field from the induced magnetization in A/m converted to nT
    Jr = magrem*mu_0*1e9; % Magnetic field from the remnant magnetization in A/m converted to nT
    % [magnetization ⋅ mu_0= (A⋅m−1 ⋅ m⋅kg⋅s−2⋅A−2) = kg⋅s−2⋅A−1 = T]; 
    % 1 nT= T*10^-9.  Multiply by 10^9 to convert T to nT

    %% Run Codes With Polygon
    % T = total field (nT)
    % I = total inclination (deg)
    % D = total declination (deg)
    % J = total magnetization (induced and remnant) in nT
[T,I,D,J] = Magnetic_Modelling_2D_Equations(Polygon,x,mu_0,Jr,Ir,Dr,Ji,Ii,Di,C);

%% magnetic field direction
x_m = 51.7;
y_m = 31;
dx = 7;
dy = 4;
x1 = [x_m - dx*cosd(I) , x_m + dx*cosd(I)];
y = [y_m + dy*sind(I) , y_m - dy*sind(I)];
dim = [0.6,0.3,0.2,0.1];
str = {'Total magnetization','J = ' + string(J) + ' (nT)','I  = '+string(I),'D = ' string(D)};

%%
figure (1)

subplot(2,1,1)

plot(x,T,'-p')
title('Magnetic Modelling KHLA')
xlim([0 100]); ylim([(min(T)-abs(min(T)*0.5)) (max(T)+abs(max(T)*0.2))]);
ylabel('Total Magnetic Intensity (nT)')
grid on

subplot(2,1,2)

plot([Polygon(:,1); Polygon(1,1)],[Polygon(:,2); Polygon(1,2)]);
set(gca,'Ydir', 'reverse'); 
axis([0 100 0 50]);
title('Polygon')
xlabel('Distance along profile (m)');
ylabel('Depth (m)')
grid on
% annotation('textarrow',x1/100,y/100); % Illustration of inclination (only when D=0 deg)
annotation('textbox',dim,'String',str,'FitBoxToText','on');
