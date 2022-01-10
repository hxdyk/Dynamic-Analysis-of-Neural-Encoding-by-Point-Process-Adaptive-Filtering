function [t_back_forth, theta,ISI] = Function_PassByPass( t_spike, T, func_x, lambda, distance, speed )

T_pass = 2*distance/speed;

t_back_forth = 0:T_pass:T;

theta = zeros(3,length(t_back_forth ));

distance_bin = 1;
x_b = distance_bin:distance_bin:distance;
for i = 2:length(t_back_forth)
    N_b = zeros(size(x_b));
    t_spike_round = t_spike(and(t_spike<=t_back_forth(i),t_spike>t_back_forth(i-1)));
    for t = t_spike_round
       N_b(round(func_x(t))+1) = N_b(round(func_x(t))+1) + 1; 
    end
    if(sum(N_b)==0)
        theta(1,i) = 0;
        theta(2,i) = 0;
        theta(3,i) = 0;
        continue;
    end
    theta(1,i) = log(sum(N_b)/T_pass);
    theta(2,i) = sum(x_b.*N_b)/sum(N_b);
    theta(3,i) = sqrt(sum((x_b-theta(2,i)).^2.*N_b)/sum(N_b));
    disp(strcat('Pass-by-Pass iteration: ',num2str(round((i-1)/(length(t_back_forth)-1)*100)) ,'%'));
end


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

theta_stat = zeros(3,length(t_observation));
for i = 1:1:length(t_back_forth)-1
    index = find(and( t_observation >= t_back_forth(i), ...
                      t_observation < t_back_forth(i+1) ) );
    theta_stat(:,index) = theta(:,i)*ones(1,length(index));
end
index = find(t_observation >= t_back_forth(end));
theta_stat(:,index) = theta(:,i)*ones(1,length(index));

rec_lambda = zeros(1,length(dN));
for i = 1:1:length(dN)
   rec_lambda(i) = lambda(x(i), theta_stat(:,i)); 
end

i_spike = find(dN==1);
ISI = zeros(length(i_spike)-1,1);
for i = 1:length(ISI)
    ISI(i) = sum(rec_lambda(i_spike(i):i_spike(i+1)))*dt_observation;
end
% typos in original text
ISI = 1-exp(-ISI);

end

