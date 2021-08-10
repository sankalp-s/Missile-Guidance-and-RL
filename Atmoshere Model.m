%clear all; clc; close all;
V = 1000;   
gamma = 60;  
g = 9.81; 
psi = 0;  
ax = 0;    ay =0 ;    az=0;  
ts = 0;    tf = 100;
syms t
d_gamma=(-(az + g*cosd(gamma)))/V;
gamma = int(d_gamma,t);
d_psi = @(t) ay/(V*cosd(gamma));
%gamma = int(F,t); %psi =  int(d_psi,t);        
Vx = @(t) V*sind(gamma);
Vy = @(t) V*sind(psi)*cosd(gamma);
Vz = @(t) V*sind(gamma)*cosd(psi);
a = @(t) ax - g*sind(gamma);
r = 0.15; S = 4*pi*(r^3)/3; cd = 0.6; m=100;


%X = RK4(Vx,0,10,0,1);
    
%% WHILE CONSIDERING ATMOSHERE

% f=input('enter the final height : ');
hf = 50000;
syms T P h rho
i=1;
plotdata(i,1)=h;
plotdata(i,2)=T;
plotdata(i,3)=P;
plotdata(i,4)=rho;
%i=i+1;
%h=0:hf;
for h=0:100:hf
    if h<=11000
        T =   15.04 - 0.00649*h;
        P=    101.29*((T+273.1)/288.08)^5.256;
        rho = P/(1718*(T+459.7));
    end
    if (11000<h)&&(h<=25000)
        T=-56.46;
        P = 22.65*exp(1.73-0.000157*h);
        rho = P/(1718*(T+459.7));
    end
    if h>25000
        T =   -131.21 + 0.00299*h;
        P=    2.488*((T+273.1)/216.6)^11.388;
        rho = P/(1718*(T+459.7));
    end
   plotdata(i,1)=h;
   plotdata(i,2)=T;
   plotdata(i,3)=P;
   plotdata(i,4)=rho;
    i=i+1;
   
end

figure(1); plot(plotdata(:,2),plotdata(:,1)); grid on; hold on;
figure(2); plot(plotdata(:,3),plotdata(:,1)); grid on;hold on;