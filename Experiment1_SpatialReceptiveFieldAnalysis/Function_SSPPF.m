function [t_observation, theta, W, ISI ] = Function_SSPPF(t_spike, T, func_x, lambda, grad, hessan)
%FUNCTION_SSPPF 此处显示有关此函数的摘要
%   此处显示详细说明

F = eye(3);
Q = diag([10^(-5), 10^(-3), 10^(-4)]);
theta_0 = [2.3;250;3.46];
W_0 = 0.001*eye(3);

% observation config
dt_observation = 0.02;
t_observation = 0:dt_observation :T;
dN = zeros(size(t_observation));
for t = t_spike
%     dN(round(t/dt_observation)+1) = dN(round(t/dt_observation)+1) + 1;
    dN(round(t/dt_observation)+1) = 1;
end
% N = cumsum(dN);

x = func_x(t_observation);

theta = zeros(length(theta_0),length(dN));
W = zeros(size(W_0,1),size(W_0,2),length(dN));
theta(:,1) = theta_0;
W(:,:,1) = W_0; % how to choose the initial value?

rec_lambda = zeros(1,length(dN));
% SSPPF
for i =2:1:length(dN)
    theta_est = F*theta(:,i-1);
    W_est = F * W(:,:,i-1) * F' + Q;
    W(:,:,i) = ( W_est^(-1) + ( grad(x(i),theta_est)'*...
                                lambda(x(i), theta_est)*...
                                dt_observation*...
                                grad(x(i),theta_est)...
                                - ...
                                (dN(i)-lambda(x(i), theta_est)*dt_observation)*...
                                hessan(x(i),theta_est)...
                               )...
                )^(-1);
    theta(:,i) = theta_est + W(:,:,i)*(  grad(x(i),theta_est)'*...
                                         (dN(i) - lambda(x(i), theta_est)*dt_observation)...
                                      );
    disp(strcat('SSPPF iteration: ',num2str(round((i-1)/(length(t_observation)-1)*100)) ,'%'));
    rec_lambda(i) = lambda(x(i), theta(:,i));
end

i_spike = find(dN==1);
ISI = zeros(length(i_spike)-1,1);
for i = 1:length(ISI)
    ISI(i) = sum(rec_lambda(i_spike(i):i_spike(i+1)))*dt_observation;
end
% typos in original text
ISI = 1-exp(-ISI);

end

