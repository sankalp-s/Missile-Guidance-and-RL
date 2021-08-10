clear all;
syms x y
f_xy = @(x,y) sin(x.^2 + abs(y));
h = 0.01;
y = 0;
x = 0:h:10;
i = 0:h:10;
%Y = zeros(1:length(x)-1);
%Y(1)=0;
for i = 0:h:10
    k1 = h*f_xy(x,y);
    k2 = h*f_xy(x + h/2 , y + k1/2);
    k3 = h*f_xy(x + h/2 , y + k2/2);
    k4 = h*f_xy(x + h , y + k3/2);
    y =  y + k1/6 + k2/3 + k3/3 + k4/6;
    f = [x;y];
    disp(f)
   % Y(i) = y;
end    
plot(x,y)
%z = max(y);