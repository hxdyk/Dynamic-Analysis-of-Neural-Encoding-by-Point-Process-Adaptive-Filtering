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

func_lambda =  @(x, alpha, mu, sigma) exp(alpha - ((x-mu).^2)./(2.*sigma.^2));
% func_lambda =  @(params) exp(params(2) - ((params(1)-params(3)).^2)./(2.*params(4).^2));

% spike data generation
J = 10000000;
dt = T/J;
t_simulation = dt:dt:T;
% time index: 1, ... , J
x = func_x(t_simulation);

% linear evolve 
alpha_linear = func_alpha_linear(t_simulation);
mu_linear = func_mu_linear(t_simulation);
sigma_linear = func_sigma_linear(t_simulation);
lambda_linear = func_lambda(x, alpha_linear, mu_linear, sigma_linear);

flag = 1;
sum = 0;
index_spike_linear = [];
threshold = 0.01;
i = 1;
exit_flag = 0;
while 1
    sum = 0;
    tao = exprnd(1);
    while abs(sum-tao) >= threshold
        sum = sum + lambda_linear(i)*dt;
        i = i+1;
        if i>J
            exit_flag = 1;
            break;
        end
    end
    if exit_flag
       break; 
    end
    index_spike_linear = [index_spike_linear, i];
end
t_spike_linear = t_simulation(index_spike_linear);

figure(1);
plot(t_spike_linear, x(index_spike_linear),'.');
hold on;
plot(t_simulation,x,'-');


% jump evolve 
alpha_jump = func_alpha_jump(t_simulation);
mu_jump = func_mu_jump(t_simulation);
sigma_jump = func_sigma_jump(t_simulation);
lambda_jump = func_lambda(x, alpha_jump, mu_jump, sigma_jump);

flag = 1;
sum = 0;
index_spike_jump = [];
threshold = 0.01;
i = 1;
exit_flag = 0;
while 1
    sum = 0;
    tao = exprnd(1);
    while abs(sum-tao) >= threshold
        sum = sum + lambda_jump(i)*dt;
        i = i+1;
        if i>J
            exit_flag = 1;
            break;
        end
    end
    if exit_flag
       break; 
    end
    index_spike_jump = [index_spike_jump, i];
end
t_spike_jump = t_simulation(index_spike_jump);
figure(2);
plot(t_spike_jump, x(index_spike_jump),'.');
hold on;
plot(t_simulation,x,'-');

% save simulation data
func_lambda_proto = func_lambda;
func_lambda = @(x,theta) func_lambda_proto(x,theta(1),theta(2),theta(3));
save('dataSpikeSimulation.mat','speed', ... 
                               'distance', ...
                               'T', ...
                               'func_x', ...
                               'func_lambda',...
                               't_spike_linear','t_spike_jump', ...
                               'J', ...
                               'func_alpha_linear', ...
                               'func_mu_linear', ...
                               'func_sigma_linear', ...
                               'func_alpha_jump', ...
                               'func_mu_jump', ...
                               'func_sigma_jump');