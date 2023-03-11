% This program is used for the article "Model-Based Design for Reactors of 
% the Modular Multilevel Converter" IEEE Transactions on Power Electronics
%
% Calculate the feasible region of Leqac-Leqdc considering short circuit
%
% Authors: Yi Zhang, Yi Xu, Maryam Saeedifard
% Contact author: Yi Zhang, yiz@ieee.org
% Date: 10-13-2022
% version: v0

clear;
clc;
close all

%% Input Parameters
% Ratings of the MMC
Udc = 60; %kV
Ug = 28.3; %kV
Ig = 1.41;%kA
Idc = 1;%kA
I0 = 1/3*Idc+1/2*Ig;

% Protection speeds
t1 = 1.07; %ms
t2 = 50; %ms

% Device short circuit capabilities
I2t = 405; %k A^2 s for diode
Isc = 5.2; %kA for IGBT

%Frequency
f = 50;%Hz
w = 2*pi*f;

%% Equations to calculate the feasible region of Leqac-Leqdc
x = 0.012:0.0001:0.2;

a1 = 27/(3*t2/1000+t1/1000);
a2 = 27*t1/1000*(4*t2/1000+t1/1000)/(4*(3*t2/1000+t1/1000)^2);
a3 = 9*(2*t2/1000+t1/1000)/(2*(3*t2/1000+t1/1000));
m = (w*t2/1000-sin(w*t2/1000))/(6*w*t2/1000+sin(2*w*t2/1000)-8*sin(w*t2/1000));
n = 1/(6*w*t2/1000+sin(2*w*t2/1000)-8*sin(w*t2/1000));
q = Udc*t1/1000./x;
b1 = 64/9*m^2-16*w*(3*t2/1000+t1/1000)*n/27;
b2 = 128/3*m^2-16*w*(2*t2/1000+t1/1000)*n/3;
b3 = 64*m^2-16*w*(t2/1000+t1/1000)*n;
b4 = 16*n;
b5 = 8*m;
Ystart = 0;
Yend = 100;

L1 = Udc*t1/1000/(3*(Isc-I0));
L2 = Udc*t1/1000/(sqrt(a1*I2t/1000-a2*I0^2)-a3*I0);
Leqac = Ug./(w*(sqrt(b1*q.^2+b2*I0.*q+b3*I0^2+b4*w*I2t/1000)-b5*(I0+q/3)));

plot([L1,L1]*1000,[Ystart Yend],'color',[0 0.498 0]);
hold on
plot([L2,L2]*1000,[Ystart Yend],'color',[0.851 0.325 0.098]);
hold on
plot(x*1000,Leqac*1000,'color',[0 0.45 0.74]);

ylim([Ystart Yend]);

title('Feasible range of the limiting reactors ');
xlabel('Leqdc (mH)');
ylabel('Leqac (mH)');

legend({'$\it{L}\rm{_{eqdc1}}\ge\frac{\it{U}\rm{_{dc}}}{\it{\lambda}\rm{_{dc1}}}$',...
                '$\it{L}\rm{_{eqdc2}}>\frac{\it{U}\rm{_{dc}}}{\it{\lambda}\rm{_{dc2}}}$',...
                '$\it{L}\rm{_{eqac}}\ge\frac{\it{\hat{U}}\rm{_{g}}}{\it{\lambda}\rm{_{ac}}}$'},...
                'Interpreter','latex','FontSize',10);