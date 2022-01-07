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

syms S E I A R
eqn1 =  -S*(params(3)*A/params(2)+params(4)*I/params(2)+params(9)*params(1))==0;
eqn2 = (params(3)*A+params(4)*I)*S/params(2)-params(6)*E-params(10)*params(1)*E==0;
eqn3 =  params(5)*params(6)*E-params(8)*I-params(12)*params(1)*I==0;
eqn4 = (1-params(5))*params(6)*E-params(7)*A-params(11)*params(1)*A==0;
eqn5 = params(7)*A+params(11)*params(1)*A+params(12)*params(1)*I+params(8)*I+params(9)*params(1)*S+params(10)*params(1)*E==0;
eqn6 = S+E+I+A+R==N;
eqn7 = S>=0;
eqn8 = E>=0;
eqn9 = I>=0;
eqn10 = A>=0;
eqn11 = R>=0;

sol = solve([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10,eqn11], [S,E,I,A,R]);
Sol_S = sol.S;
Sol_E = sol.E;
Sol_I = sol.I;
Sol_A = sol.A;
Sol_R = sol.R;


