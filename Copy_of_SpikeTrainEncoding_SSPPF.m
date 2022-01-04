clc;clear;close all;

load('dataSpikeSimulation.mat');

lambda = func_lambda;
grad = @(x,theta) [1, theta(3)^(-2)*(x-theta(2)), theta(3)^(-3)*(x-theta(2))^2]';
hessan = @(x,theta) [ 0  0                               0; ...
                      0  -theta(3)^(-2)                  -2*theta(3)^(-3)*(x-theta(2)); ...
                      0  -2*theta(3)^(-3)*(x-theta(2))   -3*theta(3)^(-4)*(x-theta(2))^2];

% observation config
dt_observation = 0.02;
t_observation = 0:dt_observation :T;
x = func_x(t_observation);

% data preprocessing
t_spike = {t_spike_linear};

dN = zeros(length(t_spike),length(t_observation));
N = zeros(length(t_spike),length(t_observation));

for  i = 1:length(t_spike)
    for t = t_spike{i}
        dN(i,round(t/dt_observation)+1) = 1;
    end
    N(i,:) = cumsum(dN(i,:));
end

% SSSPF
F = eye(3);
Q = diag([10^(-5), 10^(-3), 10^(-4)]);
theta_0 = [2.3;250;3.46];
W_0 = 0.1*eye(3);
dN = dN([1],:);

theta = zeros(length(theta_0),length(t_observation));
W = zeros(size(W_0,1),size(W_0,2),length(t_observation));
theta(:,1) = theta_0;
W(:,:,1) = W_0; % how to choose the initial value?
rec_l = zeros(4,length(t_observation));

for i =2:1:length(t_observation)
    theta_est = F*theta(:,i-1);
    W_est = F * W(:,:,i-1) * F' + Q;
    
    temp_W = zeros(size(W_0));
    err = dN(:,i)-lambda(x(i),theta_est)*dt_observation;
    temp_W = temp_W + err* hessan(x(i),theta_est);
    
    W(:,:,i) = ( W_est^(-1) + grad(x(i),theta_est)*...
                              diag(lambda(x(i),theta_est)*dt_observation)*...
                              grad(x(i),theta_est)'...
                              - temp_W ...
                )^(-1);
    theta(:,i) = theta_est + W(:,:,i)*( grad(x(i),theta_est)*err );
    rec_l(:,i) = lambda(x(i),theta_est);
    i
end

figure(3);
subplot(2,3,1);
plot(t_observation,theta(1,:));
subplot(2,3,2);
plot(t_observation,theta(2,:));
subplot(2,3,3);
plot(t_observation,theta(3,:));
subplot(2,3,4);
plot(t_observation,theta(4,:));
subplot(2,3,5);
plot(t_observation,theta(5,:));

% dN = zeros(size(t_observation));
% for t = t_spike
% %     dN(round(t/dt_observation)+1) = dN(round(t/dt_observation)+1) + 1;
%     dN(round(t/dt_observation)+1) = 1;
% end
% N = cumsum(dN);
% 
% % SSSPF
% F = eye(3);
% Q = diag([10^(-5), 10^(-3), 10^(-4)]);
% theta_0 = [2.3;250;3.46];
% W_0 = 0.1*eye(3);
% 
% theta = zeros(length(theta_0),length(t_observation));
% W = zeros(size(W_0,1),size(W_0,2),length(t_observation));
% theta(:,1) = theta_0;
% W(:,:,1) = W_0; % how to choose the initial value?
% 
% % SSPPF
% for i =2:1:length(t_observation)
%     theta_est = F*theta(:,i-1);
%     W_est = F * W(:,:,i-1) * F' + Q;
%     W(:,:,i) = ( W_est^(-1) + ( grad(x(i),theta_est)'*...
%                                 lambda(x(i),theta_est)*...
%                                 dt_observation*...
%                                 grad(x(i),theta_est)...
%                                 - ...
%                                 (dN(i)-lambda(x(i),theta_est)*dt_observation)*...
%                                 hessan(x(i),theta_est)...
%                                )...
%                 )^(-1);
%     theta(:,i) = theta_est + W(:,:,i)*(  grad(x(i),theta_est)'*...
%                                          (dN(i) - lambda(x(i),theta_est)*dt_observation)...
%                                       );
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