function  [t_observation, theta,ISI ] = Function_SDPPF(t_spike, T, func_x, lambda, grad)

% observation config
dt_observation = 0.02;
t_observation = 0:dt_observation :T;
x = func_x(t_observation);

% data preprocessing
dN = zeros(size(t_observation));
for t = t_spike
    % there might more than one spikes in an observation interval.
    dN(round(t/dt_observation)+1) = 1; 
end

% SDPPF
F = eye(3);
Q = diag([10^(-5), 10^(-3), 10^(-4)]);
theta_0 = [2.3;250;3.46];

theta = zeros(length(theta_0),length(t_observation));
theta(:,1) = theta_0;
epsilon = diag([0.02 10 1]);

rec_lambda = zeros(1,length(dN));
for i =2:1:length(t_observation)
    theta(:,i) = theta(:,i-1) + epsilon * grad(x(i),theta(:,i-1))'* (dN(i)-lambda(x(i),theta(:,i-1))*dt_observation);
    disp(strcat('SDPPF iteration: ',num2str(round((i-1)/(length(t_observation)-1)*100)) ,'%'));
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
