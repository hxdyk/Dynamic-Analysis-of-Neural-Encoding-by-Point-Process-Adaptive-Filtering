clc;clear;close all;

% init params set
speed = 125; % unit: cm/s
distance = 300; % unit: cm

T = 800; % t in (0, T]. unit: s

alpha_start = log(10);
alpha_final = log(30);
mu_start = 250;
mu_final = 150;
sigma_start = 12^(0.5);
sigma_final = 20^(0.5);

% anonymous function
func_x = @(t) (rem(speed.*t,2*distance)>distance).* ...
                                (2*distance-rem(speed.*t,2*distance)) ...
              + (rem(speed.*t,2*distance)<=distance).* ...
                                (rem(speed.*t,2*distance)) ; 
func_linear = @(t,s,f) s + (f-s).*t/T;
func_alpha_linear = @(t) func_linear(t,alpha_start,alpha_final);
func_mu_linear = @(t) func_linear(t,mu_start,mu_final);
func_sigma_linear = @(t) func_linear(t,sigma_start,sigma_final);

func_jump = @(t,s,f) (t<T/2).*s + (t>=T/2).*f;
func_alpha_jump = @(t) func_jump(t,alpha_start,alpha_final);
func_mu_jump = @(t) func_jump(t,mu_start,mu_final);
func_sigma_jump = @(t) func_jump(t,sigma_start,sigma_final);

func_lambda =  @(x,theta) exp(theta(1) - ((x-theta(2)).^2)./(2.*theta(3).^2));


% spike data generation
J = 100000;
delta_t_simulation = T/J;
t_simulation = delta_t_simulation:delta_t_simulation:T;
% time index: 1, ... , J
x = func_x(t_simulation);

% linear evolve 
alpha_linear = func_alpha_linear(t_simulation);
mu_linear = func_mu_linear(t_simulation);
sigma_linear = func_sigma_linear(t_simulation);
theta_linear = [alpha_linear;mu_linear;sigma_linear];

u = zeros(1,J);
for i = 1:1:J
    lambda = func_lambda( x(i), theta_linear(:,i));
    u(i) = binornd(1, lambda*delta_t_simulation);
end
index_spike_linear = find(u~=0);
t_spike_linear = t_simulation(index_spike_linear);
figure(1);
plot(t_spike_linear, x(index_spike_linear),'.');
hold on;
plot(t_simulation,x,'-');


% jump evolve 
alpha_jump = func_alpha_jump(t_simulation);
mu_jump = func_mu_jump(t_simulation);
sigma_jump = func_sigma_jump(t_simulation);
theta_jump = [alpha_jump;mu_jump;sigma_jump];

u = zeros(1,J);
for i = 1:1:J
    lambda = func_lambda( x(i), theta_jump(:,i));
    u(i) = binornd(1, lambda*delta_t_simulation);
end
index_spike_jump = find(u~=0);
t_spike_jump = t_simulation(index_spike_jump);
figure(2);
plot(t_spike_jump, x(index_spike_jump),'.');
hold on;
plot(t_simulation,x,'-');

% save simulation data
save('dataSpikeSimulation.mat','speed', ... 
                               'distance', ...
                               'T', ...
                               'func_x', ...
                               't_spike_linear','t_spike_jump');