clc;clear;close all;

load('dataSpikeTrainDecoding.mat');

lambda = func_lambda;
grad = @(theta) [theta(2:end)'; theta(1)*eye(length(theta)-1)];
hessan = @(theta,i) [zeros(1,i),1,zeros(1,length(theta)-i-1);[zeros(i-1,1);1;zeros(length(theta)-i-1,1)],zeros(length(theta)-1,length(theta)-1)];

% observation config
dt_observation = 0.02;
t_observation = 0:dt_observation :T;
v = func_v(t_observation);

% data preprocessing
t_spike = {t_spike_1, t_spike_2, t_spike_3, t_spike_4, t_spike_5, t_spike_6};

dN = zeros(length(t_spike),length(t_observation));
N = zeros(length(t_spike),length(t_observation));

for  i = 1:length(t_spike)
    for t = t_spike{i}
        dN(i,round(t/dt_observation)+1) = 1;
    end
    N(i,:) = cumsum(dN(i,:));
end

% SSSPF
F = diag([0.99, 1]);
Q = 10^(-5)*diag([2.5, 1]);
theta_0 = [0;3];
W_0 = zeros(2);
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
    err = dN(:,i)-lambda(theta_est)*dt_observation;
    for j = 1:length(theta_0)-1
        temp_W = temp_W + err(j)* hessan(theta_est,j);
    end
    
    W(:,:,i) = ( W_est^(-1) + grad(theta_est)*...
                              diag(lambda(theta_est)*dt_observation)*...
                              grad(theta_est)'...
                              - temp_W ...
                )^(-1);
    theta(:,i) = theta_est + W(:,:,i)*( grad(theta_est)*err );
    rec_l(:,i) = lambda(theta_est);
    i
end

figure(3);
subplot(1,2,1);
plot(t_observation,theta(1,:));
subplot(1,2,2);
plot(t_observation,theta(2,:));

% figure(3);
% subplot(2,3,1);
% plot(t_observation,theta(1,:));
% subplot(2,3,2);
% plot(t_observation,theta(2,:));
% subplot(2,3,3);
% plot(t_observation,theta(3,:));
% subplot(2,3,4);
% plot(t_observation,theta(4,:));
% subplot(2,3,5);
% plot(t_observation,theta(5,:));

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