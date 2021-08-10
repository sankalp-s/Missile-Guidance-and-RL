clear all;clc;
%close all;
syms t x
f = @(t,x) sin(t^2 + abs(x));
h = 0.001;
x = 0;
t = 0;
i = 1;
plotdata(i,1)=t;
plotdata(i,2)=x;
i = i+1;
%Y = zeros(1:length(x)-1);
%Y(1)=0;
while t < 25
    k1 = h*f(t,x);
    k2 = h*f(t + h/2 , x + k1/2);
    k3 = h*f(t + h/2 , x + k2/2);
    k4 = h*f(t + h , x + k3);
    x =  x + k1/6 + k2/3 + k3/3 + k4/6;
    
    plotdata(i,1)=t;
    plotdata(i,2)=x;
    
    t = t+h;
    i = i+1;
end

figure(1); plot(plotdata(:,1),plotdata(:,2));  hold on; grid on;
