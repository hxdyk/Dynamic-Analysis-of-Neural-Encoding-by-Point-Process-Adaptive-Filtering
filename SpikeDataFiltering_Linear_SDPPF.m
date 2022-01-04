clc;clear;close all;

load('dataSpikeSimulation.mat');

lambda = func_lambda;
grad = @(x,theta) [1, theta(3)^(-2)*(x-theta(2)), theta(3)^(-3)*(x-theta(2))^2];
hessan = @(x,theta) [ 0  0                               0; ...
                      0  -theta(3)^(-2)                  -2*theta(3)^(-3)*(x-theta(2)); ...
                      0  -2*theta(3)^(-3)*(x-theta(2))   -3*theta(3)^(-4)*(x-theta(2))^2];

% observation config
dt_observation = 0.02;
t_observation = 0:dt_observation :T;
x = func_x(t_observation);

% data preprocessing
t_spike = t_spike_jump;
dN = zeros(size(t_observation));
for t = t_spike
    % there might more than one spikes in an observation interval.
    dN(round(t/dt_observation)+1) = 1; 
end
N = cumsum(dN);

T_pass = 2*distance/speed;

t_back_forth = 0:T_pass:T;

theta = zeros(3,length(t_back_forth ));


distance_bin = 1;
x_b = distance_bin:distance_bin:distance;
N_b = zeros(size(x_b));
for i = 2:length(t_back_forth)
    N_b = zeros(size(x_b));
    t_spike_round = t_spike(and(t_spike<=t_back_forth(i),t_spike>t_back_forth(i-1)));
    for t = t_spike_round
       N_b(round(func_x(t))+1) = N_b(round(func_x(t))+1) + 1; 
    end
    theta(1,i) = log(sum(N_b)/T_pass);
    theta(2,i) = sum(x_b.*N_b)/sum(N_b);
    theta(3,i) = sqrt(sum((x_b-theta(2,i)).^2.*N_b)/sum(N_b));
end

figure(3);
subplot(1,3,1);
plot(t_back_forth,theta(1,:));
subplot(1,3,2);
plot(t_back_forth,theta(2,:));
subplot(1,3,3);
plot(t_back_forth,theta(3,:).^2);


% 
% % SDPPF
% F = eye(3);
% Q = diag([10^(-5), 10^(-3), 10^(-4)]);
% theta_0 = [2.3;250;3.46];
% 
% theta = zeros(length(theta_0),length(t_observation));
% theta(:,1) = theta_0;
% epsilon = diag([0.02 10 1]);
% for i =2:1:length(t_observation)
%     theta(:,i) = theta(:,i-1) + epsilon * grad(x(i),theta(:,i-1))'* (dN(i)-lambda(x(i),theta(:,i-1))*dt_observation);
%     i
% end
% 
% figure(3);
% subplot(1,3,1);
% plot(t_observation,theta(1,:));
% subplot(1,3,2);
% plot(t_observation,theta(2,:));
% subplot(1,3,3);
% plot(t_observation,theta(3,:).^2);