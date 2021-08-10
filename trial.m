clear all; clc; close all;
syms ay az

% Given 
g = 9.81; t =0; tf = 5;
m = 240;
T = 24000;
cd = 0.5;
S = 0.0707;
cl_alpha = 3.12;
v = 150; gamma = 0; psi = 0; x = 5000; y = 0; z = 0; 

[temp, a, P, rho] = atmosisa(x); % [T, A, P, RHO] = ATMOSISA(H) implements the mathematical 
%representation of the International Standard Atmosphere values for ambient temperature, 
%   pressure, density, and speed of sound for the input geopotential altitude.

% states(:,1) = [v; psi; gamma; x; y; z]; 

states = [v; psi; gamma; x; y; z];

dH_dU_int = 1000;
h = 0.1*pi/180; %defining the interval 
t = 0:h:tf; %range from 0 to 5 with inteval of h
gamma = 0:h:60*pi/180;

ay = 1*ones(length(t),1); %forming a matrix %lenght of t is 2865
az = 1*ones(length(t),1);
U = [ay,az];

iteration_number = 0;
error = 0.01;

display(t(1))

for i=1:length(t)-1 %i takes  value from 1 to 2865 
    
        states(:,i+1) = RK4_states(@state_dot ,t(i), states(:,i) , U(i,:) , h); %states 6x2865
                
end
vf = states(1,end);      %values of the last column of the state 1.
psif = states(2,end);    %values of the last column of the state 2.
gammaf = states(3,end);  %values of the last column of the state 3.
xf= states(4,end);       %values of the last column of the state 4.
yf= states(5,end);       %values of the last column of the state 5.
zf= states(6,end);       %values of the last column of the state 6.

lambda_vf = 0;
lambda_psif =0;
lambda_gammaf = 0;
lambda_xf= 0;
lambda_yf= 0;
lambda_zf= 0;

costates(:,length(t)) = [lambda_vf; lambda_psif; lambda_gammaf; lambda_xf; lambda_yf; lambda_zf];

for i=length(t):-1:2
        costates(:,i-1) = RK4_costates(@costate_dot, t(i) ,costates(:,i) , states(:, i),U(i,:), -h);
        
end
   
%-------------------------------------------------------
dH_dU_int = 0;
tau_controlupdate = 1e-8;

for i=1:length(t)-1 

m=240; g = 9.81; S = 0.0707; Cl = 3.12 ; Cd = 0.5; T = 24000;


lambda_v = costates(1,i); %value is zero 
lambda_psi = costates(2,i);%value is zero
lambda_gamma = costates(3,i); 
lambda_x= costates(4,i);
lambda_y= costates(5,i); 
lambda_z= costates(6,i);

%M = states(1,i);
%t = states(2,i);
v = states(1,i);
psi = states(2,i);
gamma = states(3,i); 
x= states(4,i);
y= states(5,i); 
z= states(6,i);

t_i = t(i);
U_i = U(i,:);

display(U(i,1))

%[temp, a, P, rho] = atmosisa(state(4,end));

Tw = T/(m*g); % Thrust per weight 
Sw = rho*a^2*S/(2*m*g); %Surface per weight



end
%------------------------------------------

function X_dot = state_dot(tau,X,U)

m=240; g = 9.81; S = 0.0707; cl_alpha = 3.12 ; cd = 0.5; T = 24000;

Tw = T/(m*g); 
t = tau;

v = X(1);
psi = X(2);
gamma = X(3);
x = X(4);
y = X(5);
z = X(6);

ay = U(1);
az = U(2);

[temp, a, P, rho] = atmosisa(5000); 
Sw = rho*a^2*S/(2*m*g);
Q = 0.5*rho*v^2; 
D = Q*cd*S;
alpha = m*sqrt(ay^2 + az^2)/(Q*S*cl_alpha);
ax = (T*cos(alpha) - D)/m;

v_dot = ax - g*sin(gamma);
gamma_dot = (-(az + g*cos(gamma)))/v;
psi_dot = ay/(v*cos(gamma));
x_dot = v*sin(gamma);
y_dot = v*cos(gamma)*sin(psi);
z_dot = v*cos(gamma)*cos(psi); 

X_dot=[v_dot; psi_dot; gamma_dot; x_dot; y_dot; z_dot];

end

function lambda_dot = costate_dot(tau,lambda,X,U)

m=240; g = 9.81; S = 0.0707; cl_alpha = 3.12 ; cd = 0.5; T = 24000;

t = tau;

lambda_v = lambda(1);
lambda_psi = lambda(2);
lambda_gamma = lambda(3);
lambda_x = lambda(4);
lambda_y = lambda(5);
lambda_z =lambda(6);

v = X(1);
psi = X(2);
gamma = X(3);
x = X(4);
y = X(5);
z = X(6);

ay = U(1);
az = U(2);

[temp, a, P, rho] = atmosisa(5000);

Tw = T/(m*g); 
Sw = rho*a^2*S/(2*m*g);

lambda_v_dot =(ay*conj(lambda_psi))/(v^2*cos(gamma)) - (conj(lambda_gamma)*(az + g*cos(gamma)))/v^2 - cos(gamma)*conj(lambda_z)*cos(psi) - cos(gamma)*conj(lambda_y)*sin(psi) - conj(lambda_x)*sin(gamma) - (conj(lambda_v)*(4*T*m*sin((2*m*(ay^2 + az^2)^(1/2))/(S*cl_alpha*rho*v^2))*(ay^2 + az^2)^(1/2) - S^2*cd*cl_alpha*rho^2*v^4))/(S*cl_alpha*m*rho*v^3);
lambda_psi_dot = v*cos(gamma)*conj(lambda_z)*sin(psi) - v*cos(gamma)*conj(lambda_y)*cos(psi);
lambda_gamma_dot = g*cos(gamma)*conj(lambda_v) - v*cos(gamma)*conj(lambda_x) - (g*conj(lambda_gamma)*sin(gamma))/v + v*conj(lambda_z)*cos(psi)*sin(gamma) + v*conj(lambda_y)*sin(gamma)*sin(psi) - (ay*conj(lambda_psi)*sin(gamma))/(v*cos(gamma)^2);
lambda_x_dot = 0;
lambda_y_dot = 0;
lambda_z_dot = 0;

lambda_dot = [lambda_v_dot; lambda_psi_dot; lambda_gamma_dot; lambda_x_dot; lambda_y_dot; lambda_z_dot];

end

% Rk4 on states 
function y = RK4_states(f,tau,X,U,h)
   
    k1 = h*f(tau,X,U); 
    k2 = h*f(tau + 0.5*h , X + 0.5*k1 , U);
    k3 = h*f(tau + 0.5*h , X + 0.5*k2 , U);
    k4 = h*f(tau + h , X+ k3 , U);
    y = X + (k1 + 2*k2 + 2*k3 + k4)/6;

end

%RK4 on costates
function y = RK4_costates(f,tau,lambda,X,U,h)
   
    k1 = h*f(tau,lambda,X,U); 
    k2 = h*f(tau + 0.5*h , lambda + 0.5*k1 , X , U);
    k3 = h*f(tau + 0.5*h , lambda + 0.5*k2 , X  , U);
    k4 = h*f(tau + h ,lambda + k3, X , U);
    y = lambda + (k1 + 2*k2 + 2*k3 + k4)/6;

end