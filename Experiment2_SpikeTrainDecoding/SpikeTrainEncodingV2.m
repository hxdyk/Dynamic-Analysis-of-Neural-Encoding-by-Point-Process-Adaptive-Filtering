clc;clear;close all;

T = 800; % unit: s

% Generate simulation data.
func_beta_1 = @(t) 3*ones(size(t));
func_beta_2 = @(t) -3*ones(size(t));
func_beta_3 = @(t) (0<t & t<=200).*0 + ...
              (200<t & t<=400).*0.0125.*(t-200) + ...
              (400<t & t<=800).*2.5;
func_beta_4 = @(t) (0<t & t<=400).*0 - ...
              (400<t & t<=600).*0.0125.*(t-400) - ...
              (600<t & t<=800).*2.5;
func_beta_5 = @(t) 2.5*ones(size(t));
func_beta_6 = @(t) -2.5*ones(size(t));

exp_mu = 1;

func_lambda_proto = @(beta,v) exp_mu*exp(beta.*v);
func_lambda = @(theta) exp_mu*exp(theta(2:end)*theta(1));

t_step = 0.001; % Unit: s
t_simulation = 0:t_step:T;

var_ita_v = 2.5e-5;

dv_sim = normrnd(0, sqrt(var_ita_v),size(t_simulation));
v_sim = cumsum(dv_sim);

v_sim = v_sim/max(abs(v_sim));

func_v = @(t) (0<=t & t<=T).*v_sim(round(t./t_step)+1); % t with unit s, look up for real time velocity.

J = T/t_step;
dt = T/J;
t_simulation = dt:dt:T;
% time index: 1, ... , J
v = v_sim(1,2:end);
beta_1 = func_beta_1(t_simulation);
beta_2 = func_beta_2(t_simulation);
beta_3 = func_beta_3(t_simulation);
beta_4 = func_beta_4(t_simulation);
beta_5 = func_beta_5(t_simulation);
beta_6 = func_beta_6(t_simulation);

lambda_1 = func_lambda_proto(beta_1,v);
lambda_2 = func_lambda_proto(beta_2,v);
lambda_3 = func_lambda_proto(beta_3,v);
lambda_4 = func_lambda_proto(beta_4,v);
lambda_5 = func_lambda_proto(beta_5,v);
lambda_6 = func_lambda_proto(beta_6,v);

lambda = {lambda_1,lambda_2,lambda_3,lambda_4,lambda_5,lambda_6};
t_spike = cell(1,6);

for k = 1:6
    flag = 1;
    sum = 0;
    index_spike = [];
    threshold = 0.01;
    i = 1;
    exit_flag = 0;
    while 1
        sum = 0;
        tao = exprnd(1);
        while abs(sum-tao) >= threshold
            sum = sum + lambda{1,k}(i)*dt;
            i = i+1;
            if i>J
                exit_flag = 1;
                break;
            end
        end
        if exit_flag
           break; 
        end
        index_spike = [index_spike, i];
    end
    t_spike{1,k} = t_simulation(index_spike);
    
end

t_spike_1 = t_spike{1,1};
t_spike_2 = t_spike{1,2};
t_spike_3 = t_spike{1,3};
t_spike_4 = t_spike{1,4};
t_spike_5 = t_spike{1,5};
t_spike_6 = t_spike{1,6};

% save simulation data
save('dataSpikeTrainDecoding.mat','func_lambda',...
                                  'func_beta_1', ...
                                  'func_beta_2', ...
                                  'func_beta_3', ...
                                  'func_beta_4', ...
                                  'func_beta_5', ...
                                  'func_beta_6', ...
                                  't_spike_1', ...
                                  't_spike_2', ...
                                  't_spike_3', ...
                                  't_spike_4', ...
                                  't_spike_5', ...
                                  't_spike_6', ...
                                  'func_v',...
                                  'v_sim',...
                                  'T');
