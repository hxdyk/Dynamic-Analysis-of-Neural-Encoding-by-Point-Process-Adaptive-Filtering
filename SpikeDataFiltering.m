clc;clear;close all;

load('dataSpikeSimulation.mat');

% observation config
delta_t_observation = 0.02;
t_observation = delta_t_observation :delta_t_observation :T;
x = func_x(t_observation);

% data preprocessing
t_spike = t_spike_linear;
dN = zeros(size(t_observation));
for i = 1:1:numel(t_spike)
    dN(ceil(t_spike(i))) = dN(ceil(t_spike(i))) + 1;
end
N = cumsum(dN);

% SSSPF
F = eye(3);
Q = diag([10^(-5), 10^(-3), 10^(-4)]);
theta = [0;0;0];
for i =1:1:size(t_observation)
    theta_est = 
end
