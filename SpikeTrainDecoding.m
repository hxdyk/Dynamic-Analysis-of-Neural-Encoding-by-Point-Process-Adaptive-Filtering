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

mu = log(1);

func_lambda = @(theta) exp(mu+theta(2:end)*theta(1));

t_step = 0.001; % Unit: s
t_simulation = 0:t_step:T;

var_ita_v = 2.5e-5;
v_sim = zeros(size(t_simulation));
for i = 2:length(t_simulation)
    v_sim(i) = v_sim(i-1) + normrnd(0, var_ita_v);
end
func_v = @(t) (0<=t & t<=T).*v_sim(round(t./t_step)+1); % t with unit s, look up for real time velocity.

J = 40000;
delta_t_simulation = T/J;
t_simulation = delta_t_simulation:delta_t_simulation:T;
% time index: 1, ... , J
v = func_v(t_simulation);
beta_1 = func_beta_1(t_simulation);
beta_2 = func_beta_2(t_simulation);
beta_3 = func_beta_3(t_simulation);
beta_4 = func_beta_4(t_simulation);
beta_5 = func_beta_5(t_simulation);
beta_6 = func_beta_6(t_simulation);

u = zeros(6,length(t_simulation));
rec_lambda = zeros(6,length(t_simulation));
for i = 1:1:J
    lambda = func_lambda([v(i) beta_1(i) beta_2(i) beta_3(i) beta_4(i) beta_5(i) beta_6(i)]);
    u(:,i) = binornd(1, lambda*delta_t_simulation)';
    rec_lambda(:,i) = lambda';
end
t_spike_1 = t_simulation(u(1,:)~=0);
t_spike_2 = t_simulation(u(2,:)~=0);
t_spike_3 = t_simulation(u(3,:)~=0);
t_spike_4 = t_simulation(u(4,:)~=0);
t_spike_5 = t_simulation(u(5,:)~=0);
t_spike_6 = t_simulation(u(6,:)~=0);

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
                                  'mu',...
                                  'T');
                            
%%

% plot(t_simulation,rec_lambda(1,:));
% hold on;plot(t_simulation,rec_lambda(2,:));
% hold on;plot(t_simulation,rec_lambda(3,:),'*-');
% hold on;plot(t_simulation,rec_lambda(4,:));
% hold on;plot(t_simulation,rec_lambda(5,:));
% hold on;plot(t_simulation,rec_lambda(6,:));
% 
% plot(t_simulation,cumsum(u(1,:)) );
% hold on;plot(t_simulation,cumsum(u(2,:)));
% hold on;plot(t_simulation,cumsum(u(3,:)));
% hold on;plot(t_simulation,cumsum(u(4,:)));
% hold on;plot(t_simulation,cumsum(u(5,:)));
% hold on;plot(t_simulation,cumsum(u(6,:)));
% legend('1','2','3','4','5','6');
% 
% plot(t_observation,N(1,:) );
% hold on;plot(t_observation,N(2,:));
% hold on;plot(t_observation,N(3,:));
% hold on;plot(t_observation,N(4,:));
% hold on;plot(t_observation,N(5,:));
% hold on;plot(t_observation,N(6,:));
% legend('1','2','3','4','5','6');