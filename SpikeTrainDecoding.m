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

func_lambda = @(theta) exp_mu*exp(theta(2:end)*theta(1));

t_step = 0.001; % Unit: s
t_simulation = 0:t_step:T;

var_ita_v = 2.5e-5;
% drift = 0;
% std = sqrt(var_ita_v);         % std = sqrt(variance)
% pd = makedist('Normal',drift,std);
% nsteps = T/t_step;
% Z = random(pd,nsteps,1);
% X = [0; cumsum(Z)];
% plot(0:nsteps,X)          % alternatively:  stairs(0:nsteps,X)   

dv_sim = normrnd(0, sqrt(var_ita_v),size(t_simulation));
v_sim = cumsum(dv_sim);
% v_sim = zeros(size(t_simulation));
% for i = 2:length(t_simulation)
%     v_sim(i) = dv_sim(i) + 0.99*v_sim(i-1);
% end
v_sim = v_sim/max(abs(v_sim));

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
for i = 1:1:length(t_simulation)
    lambda = func_lambda([v(i) beta_1(i) beta_2(i) beta_3(i) beta_4(i) beta_5(i) beta_6(i)]);
    u(:,i) = binornd(1, lambda*delta_t_simulation)';
%     u(:,i) = binornd(1, lambda*t_step)';
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
                                  'rec_lambda',...
                                  'u',...
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