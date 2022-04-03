clc;clear;close all;

index_experiment = 2; % 1 for data_linear, 2 for data_jump.
index_param = [1,2,3]; % 1 for alpha, 2 for mu, 3 for sigma.
index_method = [1,2,3,4]; % 1 for Pass-By-Pass, 2 for EKF, 3 for SDPPF, 4 for SSPPF
%% Load data
load('dataSpikeSimulation.mat');

lambda =  @(x,theta) exp(theta(1) - ((x-theta(2)).^2)./(2.*theta(3).^2));
grad = @(x,theta) [1, theta(3)^(-2)*(x-theta(2)), theta(3)^(-3)*(x-theta(2))^2];
hessan = @(x,theta) [ 0  0                               0; ...
                      0  -theta(3)^(-2)                  -2*theta(3)^(-3)*(x-theta(2)); ...
                      0  -2*theta(3)^(-3)*(x-theta(2))   -3*theta(3)^(-4)*(x-theta(2))^2];


t_spike_data = {t_spike_linear,t_spike_jump};
t_spike = t_spike_data{index_experiment};

%% Data processing
% data initialization
t_obs_PBP = []; theta_PBP = [];
t_obs_EKF = []; theta_EKF = []; W_EKF = [];
t_obs_SDPPF = []; theta_SDPPF = [];
t_obs_SSPPF = []; theta_SSPPF = []; W_SSPPF = [];

for i = index_method
    switch i
        case 1
            [t_obs_PBP, theta_PBP,ISI_PBP] = Function_PassByPass( t_spike, T, func_x, lambda, distance, speed );
        case 2
            [t_obs_EKF, theta_EKF, W_EKF, ISI_EKF ] = Function_EKF(t_spike, T, func_x, lambda, grad, hessan);
        case 3
            [t_obs_SDPPF, theta_SDPPF, ISI_SDPPF] = Function_SDPPF(t_spike, T, func_x, lambda, grad);
        case 4
            [t_obs_SSPPF, theta_SSPPF, W_SSPPF, ISI_SSPPF] = Function_SSPPF(t_spike, T, func_x, lambda, grad, hessan);
    end
end

%% Visualizing results

func_groundtruth = {func_alpha_linear, func_mu_linear, func_sigma_linear;
                    func_alpha_jump,   func_mu_jump,   func_sigma_jump};
t_obs = {t_obs_PBP, t_obs_EKF, t_obs_SDPPF, t_obs_SSPPF};
theta = {theta_PBP, theta_EKF, theta_SDPPF, theta_SSPPF};
W = {[], W_EKF, [], W_SSPPF};
title = {'Pass-by-Pass','EKF','SDPPF','SSPPF'};

figure(index_experiment);
row = length(index_param);
col = length(index_method);
MSE = zeros(row, col);
coverage = zeros(row, col);
for ii = 1:row
    i = index_param(ii);
    for j = 1:col
        subplot(row,col,j+ col*(ii-1));
        k = index_method(j);
        plot(t_obs{k}, theta{k}(i,:), t_obs{k}, func_groundtruth{index_experiment,i}(t_obs{k}));
        MSE(ii,j) = mean(( theta{k}(i,:) - func_groundtruth{index_experiment,i}(t_obs{k}) ).^2);
        if ~isempty(W{k})
            hold on;
            plot(t_obs{k}, theta{k}(i,:)+2.475*squeeze(sqrt(W{k}(i,i,:)))','b-', t_obs{k}, theta{k}(i,:)-2.475*squeeze(sqrt(W{k}(i,i,:)))','b-');
            range = and( func_groundtruth{index_experiment,i}(t_obs{k}) <= theta{k}(i,:)+2.475*squeeze(sqrt(W{k}(i,i,:)))',...
                         func_groundtruth{index_experiment,i}(t_obs{k}) >= theta{k}(i,:)-2.475*squeeze(sqrt(W{k}(i,i,:)))');
            coverage(ii,j) = sum(range)/(length(t_obs{k})-1);
        end
%         legend(title{j},'Ground Truth');
    end
end 

%%
pd_uni = makedist('Uniform');

figure(3);
ISI = {ISI_PBP,ISI_EKF;ISI_SDPPF,ISI_SSPPF};
h = cell(2,2);
KS = zeros(2,2);
switch index_experiment
    case 1
        interval = 0.044;
    case 2
        interval = 0.040;
end
for i = 1:2
    for j = 1:2
        subplot(2,2,j+(i-1)*2);
        ISI_data = ISI{i,j};
        x = linspace(0,1,length(ISI_data))';
        test_cdf = [x,cdf(pd_uni,x)];
        [h{i,j},p{i,j}] = kstest(sort(ISI_data),'CDF',test_cdf);
        KS(i,j) = max(abs(sort(ISI_data)-cdf(pd_uni,x)));
        plot(x,sort(ISI_data));
        hold on;
        plot(x,cdf(pd_uni,x),'r-',x,cdf(pd_uni,x)+interval,'r-',x,cdf(pd_uni,x)-interval,'r-');
        axis([0,1,0,1]);
    end 
end
