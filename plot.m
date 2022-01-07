%declare parameters
alpha=0.1; 
N=49520000; 
beta_1=0.8756; 
beta_2=0.7833; 
gamma=0.0348;
sigma=1/5.2;
k_1=0.2;
k_2=1/7+1/15.16;
theta_s=0.85;
theta_e=0.85;
theta_a=0.85;
theta_i=0.85;
params=[alpha, N, beta_1, beta_2, gamma, sigma, k_1, k_2, theta_s, theta_e, theta_a, theta_i];


%prepare ode for solver
odefun = @(t,y) my_ode(y,params);

t0 = 0; tf = 100;
tspan = [t0, tf]; % this tells the solver to solve from t0 to tf

%declare initial condition
S0=49519960;
E0=20;
A0=19;
I0=1;
R0=0;
Y_init=[S0; E0; I0; A0; R0];

%make solver more accurate
ODE_options = odeset('RelTol',1e-8);

%solve the ODE and store the output into the solution struct sol
sol = ode45(odefun,tspan,Y_init,ODE_options);

%plot E against t
num_columns=size(sol.y,2);
t=1:1:num_columns;
plot(t,sol.y(2,:));
xlabel('t');
ylabel('E');
title('E against t');

% %plot I aginst t
% num_columns=size(sol.y,2);
% t=1:1:num_columns;
% plot(t,sol.y(3,:));
% xlabel('t');
% ylabel('I');
% title('I against t');

function dYdt = my_ode(Y,params)
    dYdt= zeros(5,1);
    dYdt(1)= -Y(1)*(params(3)*Y(4)/params(2)+params(4)*Y(3)/params(2)+params(9)*params(1));
    dYdt(2)=(params(3)*Y(4)+params(4)*Y(3))*Y(1)/params(2)-params(6)*Y(2)-params(10)*params(1)*Y(2);
    dYdt(3)=params(5)*params(6)*Y(2)-params(8)*Y(3)-params(12)*params(1)*Y(3);
    dYdt(4)=(1-params(5))*params(6)*Y(2)-params(7)*Y(4)-params(11)*params(1)*Y(4);
    dYdt(5)=params(7)*Y(4)+params(11)*params(1)*Y(4)+params(12)*params(1)*Y(3)+params(8)*Y(3)+params(9)*params(1)*Y(1)+params(10)*params(1)*Y(2);


    
end