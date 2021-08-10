clear all; clc;
syms v psi gamma x y z ay az g m T cd rho S cl_alpha alpha %for system
syms lambda_v lambda_gamma lambda_psi lambda_x lambda_y lambda_z %for Lambda

Q = 1/2*rho*v^2; %Dynamic pressure; Q = Dynamic pressure
D = Q*S*cd; %Drag equation D = drag; cd = coefficient of Drag; S = The  wing surface.  
alpha = m*sqrt(ay^2 + az^2)/(Q*S*cl_alpha);%Alpha = angle of attack i.e angle between wing and airfoil
ax = (T*cos(alpha) - D)/m;

v_dot = ax - g*sin(gamma);
psi_dot = ay/(v*cosd(gamma));
gamma_dot = (-(az + g*cosd(gamma)))/v;
x_dot = v*sin(gamma);
y_dot = v*cos(gamma)*sin(psi);
z_dot = v*cos(gamma)*cos(psi);

%rho = a+b*exp(-(alpha*height+beta))

%Question
%Given Cost function 
% J = 1/2*(integral|0 to tf(ay^2 + az^2)) > objective function, tf = 5 seconds
% x_0 given x_f = free for v psi x y z , gamma to achieve. gamma_final = 60
% degree after 5 seconds.

% This method is called as Gradient method or Steepest Descent method.

%% Finding the State Equations (From symbolic Toolbox)
f1 = ax - g*sin(gamma) ; %v_dot
f2 = ay/(v*cos((gamma))) ; %psi_dot
f3 = -(az + g*cos((gamma)))/v; %gamma_dot
f4 = v*sin(gamma) ; %x_dot
f5 = v*cos(gamma)*sin(psi) ; %y_dot
f6 = v*cos(gamma)*cos(psi); %z_dot

%pause
% J = 1/2(int((ay^2 + az^2) , t));
 
L = 0.5*((ay^2 + az^2)); %Function in integral
f = [f1;f2;f3;f4;f5;f6]; 

lambda = [lambda_v; lambda_psi; lambda_gamma; lambda_x; lambda_y; lambda_z];

%Note: Lamba and f have to be corresponding. 

%% Forming Hamiltonian
H = L + lambda'*f;
% disp(simplify(H));


%% Forming Costate equation
lambda_v_dot = simplify(-diff(H,v));
lambda_psi_dot = simplify(-diff(H,psi));
lambda_gamma_dot = simplify(-diff(H,gamma));
lambda_x_dot = simplify(-diff(H,x));
lambda_y_dot = simplify(-diff(H,y));
lambda_z_dot = simplify(-diff(H,z));


%% Optimal Control Equations (U is the control functoin)
dH_dU_ay = simplify(diff(H,ay))
dH_dU_az = simplify(diff(H,az))
