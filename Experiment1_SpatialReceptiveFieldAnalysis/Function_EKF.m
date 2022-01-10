function [t_observation, theta, W, ISI ] = Function_EKF(t_spike, T, func_x, lambda, grad, hessan)
%FUNCTION_EKF 此处显示有关此函数的摘要
%   此处显示详细说明

% observation config
dt_observation = 0.02;
t_observation = 0:dt_observation :T;
x = func_x(t_observation);

% data preprocessing
dN = zeros(size(t_observation));
for i = 1:1:numel(t_spike)
%     dN(round(t/dt_observation)+1) = dN(round(t/dt_observation)+1) + 1;
    dN(round(t_spike(i)/dt_observation)+1) = 1;
end
% N = cumsum(dN);

width = 0.25; % unit: s
sigma_width = 0.25/3;
filter = @(x) sqrt(2/pi).*exp(-0.5*(x/sigma_width).^2)./sigma_width;

r = zeros(size(t_observation));
rec_lambda = zeros(1,length(dN));
for t = t_spike
    i_start = (round(t/dt_observation)+1);
    i_end = min((round((t+width)/dt_observation)+1), length(t_observation));
    for i_impact = i_start:i_end
        r(i_impact) = r(i_impact) + filter( t_observation(i_impact)-t );
    end
end

% axis = 0:0.001:1;
% plot(axis,filter(axis));
% plot(r,'*-');

% observation config
dt_observation = 0.02;
t_observation = 0:dt_observation :T;
x = func_x(t_observation);


F = eye(3);
Q = diag([10^(-5), 10^(-3), 10^(-4)]);
theta_0 = [2.3;250;3.46];
W_0 = 0.1*eye(3);

theta = zeros(length(theta_0),length(t_observation));
W = zeros(size(W_0,1),size(W_0,2),length(t_observation));
theta(:,1) = theta_0;
W(:,:,1) = W_0; % how to choose the initial value?

rec_lambda = zeros(1,length(dN));
for i =2:1:length(t_observation)
    theta_est = F*theta(:,i-1);
    W_est = F * W(:,:,i-1) * F' + Q;
    W(:,:,i) = ( W_est^(-1) + grad(x(i),theta_est)'*...
                              lambda(x(i),theta_est)*...
                              dt_observation*...
                              grad(x(i),theta_est) )^(-1);
    theta(:,i) = theta_est + W(:,:,i)*(  grad(x(i),theta_est)'*...
                                         (r(i) - lambda(x(i),theta_est))*...
                                         dt_observation...
                                      );
    disp(strcat('EKF iteration: ',num2str(round((i-1)/(length(t_observation)-1)*100)) ,'%'));
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

