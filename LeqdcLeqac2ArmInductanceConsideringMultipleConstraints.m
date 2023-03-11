% This program is used for the article "Model-Based Design for Reactors of 
% the Modular Multilevel Converter" IEEE Transactions on Power Electronics
%
% Consider different impacts to obtain arm inductance L0
%
% Authors: Yi Zhang, Yi Xu, Maryam Saeedifard
% Contact author: Yi Zhang, yiz@ieee.org
% Date: 10-13-2022
% version: v0

clear;
clc;
close all
%% Input parameters
% Ratings of the MMC
S = 60;%MVA
Udc = 60;%kV
Ug = 28.3;%kV
Ig = 1.41;%kA
Idc = 1;%kA
I0 = 1/3*Idc+1/2*Ig;

%Protection speeds
t1 = 1.07;%ms
t2 = 50;%ms

%Device short circuit capabilities
I2t = 405;%k A^2 s for diode
Isc = 5.2;%kA for IGBT

%Frequency
f = 50;%Hz
w = 2*pi*f;

%MMC parameters
N = 20;
Usm = Udc/N;
C0 = 2.65e-3;%F
Ls = 6.37;%Equivalent inductance of the AC grid(mH)
ma = 0.96; %Modulation index

%Harmonic limit
THDpcc = 0.015;

%Leqdc range(mH)
Xstart = 16;
Xend = 150;

%Leqac range(mH)
Ystart = 0;
Yend = 100;

%Operation
phic = pi/4;

%% Equations
% Limiting short circuit transient
x = Xstart/1000:0.0005:Xend/1000;
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

L1 = Udc*t1/1000/(3*(Isc-I0));%\nu dc1
L2 = Udc*t1/1000/(sqrt(a1*I2t/1000-a2*I0^2)-a3*I0);%\nu dc2
Leqac_dc = Ug./(w*(sqrt(b1*q.^2+b2*I0.*q+b3*I0^2+b4*w*I2t/1000)-b5*(I0+q/3)));%\nu ac
Leqac = Leqac_dc'*1000;

% Avoid current resonance
L0min = N*(3+2*ma^2)/(48*w^2*C0)*1000;

% THD
t = 0;
for k = 1:1:25
    n = 2*k+1;
    f = 0;
    for s = 1:1:min(round(ma*Udc/(2*Usm),N/2))
        theta = asin((s-0.5)*Usm*2/(ma*Udc));
        f = f+cos(n*theta);
    end
    t = t+(f/n)^2;
end
f1 =0;
for j = 1:1:min(round(ma*Udc/(2*Usm),N/2))
        thetaj = asin((j-0.5)*Usm*2/(ma*Udc));
        f1 = f1+cos(thetaj);
end

THDc = sqrt(t)/f1;
Leqacmin = THDc*Ls/THDpcc-Ls;

% Interface to AC grid
Leqacmax = 1000*9*Ug^2*(sqrt(4*S^2*w^2*Udc^2/(9*Ug^2)-16*S^2*(cos(phic))^2*w^2/9)+4*S*sin(phic)*w/3)/(8*S^2*w^2);

% Arm inductance
Leqdc = x'*1000;
i = 1;
J = ones(1);
num = length(Leqdc);
for s = 1:1:num
    y = Leqac(s,:);
    m = Leqdc(s,:);
    if y <= Leqacmax && y >= Leqacmin
        for r = m:1:Leqdc(num,:)
            L = L0min;
            for k = 1:300
                if (L <= 1.5*r) && (L <= 2*y)
                    L = L+1;
                else
                    J(i,1) = r;
                    J(i,2) = y;
                    J(i,3) = L-1;
                end
            end
            i = i+1;
        end
    else
    end
end
c = J(:,1);
d = J(:,2);
e = J(:,3);
[row,column]=size(e);
Fm = ones(row,1)*L0min;
Fxoy = ones(row,1)*0;
xtemp=linspace(min(c),max(c),100);
ytemp=linspace(min(d),max(d),100);
[X,Y]=meshgrid(xtemp,ytemp);
Top = griddata(c,d,e,X,Y,'cubic');
Flatmin = griddata(c,d,Fm,X,Y,'cubic');
Flatxoy = griddata(c,d,Fxoy,X,Y,'cubic');

% Side face
count = 1;
for g = 1:1:num
    x2 = Leqdc(g,:);
    y2 = Leqac(g,:);
    z2 = L0min;
    if y2 <= Leqacmax && y2 >= Leqacmin
        for k = 1:300
            if (z2 <= 1.5*x2) && (z2 <= 2*y2)
                z2 = z2+1;
            else
                O(count,:) = x2;
                U(count,:) = y2;
                Z(count,:) = z2-1;
            end
        end
        count = count+1;
    else
    end
end              
h = max(Z);
v = h-L0min+1;
A = repmat(O,1,round(h-L0min+1));
B = repmat(U,1,round(h-L0min+1));
data = length(O);
for n = 1:1:data
    q = L0min;
    for j = 1:1:round(h-L0min+1)
        Q(n,j) = q;
        if q < Z(n,1)
            q = q+1;
        else
        end
    end
end

%% Plot
Y1 = Ystart:Yend;
num2 = length(Y1);
X1 = ones(num2,1)*L1*1000;
X2 = ones(num2,1)*L2*1000;
P1 = zeros(num2,1);
P2 = zeros(num,1);

X3 = 0:Xend;
Ymax = ones(length(X3),1)*Leqacmax;
Ymin = ones(length(X3),1)*Leqacmin;
P3 = zeros(length(X3),1);

plot3(X1,Y1,P1,'color',[0 0.498 0]);%\nu dc1
hold on
plot3(X2,Y1,P1,'color',[0.851 0.325 0.098]);%\nu dc2
hold on
plot3(O,U,Z,'color',[0 0.45 0.74]);%\nu ac
hold on
plot3(X3,Ymax,P3,'color',[0 0.45 0.74]);%Max Leqac
hold on
plot3(X3,Ymin,P3,'color',[0 0.45 0.74]);%Min Leqac
hold on
mesh(X,Y,Top);
hold on
mesh(X,Y,Flatmin);
hold on
mesh(X,Y,Flatxoy);
hold on
mesh(A,B,Q);
colormap parula
colorbar
grid on
xlabel('\it{L}\rm{_{eqdc}}(mH)','Fontname', 'Times New Roman');
ylabel('\it{L}\rm{_{eqac}}(mH)','Fontname', 'Times New Roman');
zlabel('\it{L}\rm{_{0}}(mH)','Fontname', 'Times New Roman');

