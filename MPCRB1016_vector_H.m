%% Setup Parameters
clear; clc

%% Setup Experiment
nSamples = 50;
%% Compute PMCRB

%% MAP Loop
n_trial = 100000;

dimension = 3;
% Correctly Specified Parameters
param_ps.mu = [10;20;5];
param_ps.sigma = 0.5*eye(dimension);
param_ps.h = 1;
rho = 0.5;
Q = [1 rho rho^2;rho 1 rho;rho^2 rho 1];
param_ps.Sigma = 0.04*Q;

% Misspecified Parameters
param_th.mu = [8;18;6];
param_th.sigma = 0.5*eye(dimension);
% param_th.h = 1;
param_th.Sigma = 0.1*eye(dimension);

% Storing Results
ps = zeros(3,n_trial);
x = zeros(3,nSamples,n_trial);
est_th = zeros(3,n_trial);
%     error = zeros(n_trial,1);
theta0 = zeros(3,n_trial);

% Simulation Variable Setup
hCandidate = (0.2:0.1:3);
nCandidate = length(hCandidate);
RMSE = zeros(3,3, nCandidate);
bound = zeros(3,3, nCandidate);
% CRB = zeros(3,3, nCandidate);

for n = 1:nCandidate
    n
    param_th.h = hCandidate(n);
    parfor m = 1:n_trial
        %% Data Generation
        % Data of Prior
        ps(:,m) = mvnrnd(param_ps.mu,param_ps.sigma);
        % Data of x distribution
        x(:,:,m) = mvnrnd(param_ps.h*ps(:,m),param_ps.Sigma,nSamples)';
        % Pseudotrue
        theta0(:,m) = (nSamples*param_th.h'*param_th.Sigma^(-1)*param_th.h + param_th.sigma^(-1))^(-1)...
            *(nSamples*param_th.h'*param_th.Sigma^(-1)*param_ps.h*ps(:,m)+ param_th.sigma^(-1)*param_th.mu);
        % Estimation
        est_th(:,m) = (nSamples*param_th.h'*param_th.Sigma^(-1)*param_th.h + param_th.sigma^(-1))^(-1)...
            *(param_th.h'*param_th.Sigma^(-1)*sum(x(:,:,m),2)+ param_th.sigma^(-1)*param_th.mu); % Implement MAP
        %         mean_est_th(:,m) = mean(est_th(1,1:m)); % mean of estimate convergence
    end

    %% Error Record
    error = est_th - ps;
    covar = zeros(3,3);
    for m = 1:n_trial
        covar = error(:,m)*error(:,m)' + covar;
    end
    RMSE(:,:,n) = sqrt(covar/n_trial);

    %% Bound Computation
    A = (nSamples*param_th.h'*param_th.Sigma^(-1)*param_th.h + param_th.sigma^(-1))^(-1)...
        *(nSamples*param_th.h'*param_th.Sigma^(-1)*param_ps.h);
    B = nSamples*param_ps.h'*param_ps.Sigma^(-1)*param_ps.h + param_ps.sigma^(-1);
%     CRB(:,:,n) = A*B^(-1)*A';

    A1 = nSamples*param_th.h'*param_th.Sigma^(-1)*param_th.h + param_th.sigma^(-1);
    B1 = nSamples*param_th.h'*param_th.Sigma^(-1)*param_ps.h;
    C1 = param_th.sigma^(-1)*param_th.mu;
    D1 = param_ps.sigma + param_ps.mu*param_ps.mu';
    E1 = param_ps.mu;
    rrT = A1^(-1)*(B1*D1*B1' + B1*E1*C1' + C1*E1'*B1' + C1*C1')*(A1^(-1))'...
        - A1^(-1)*(B1*D1' + C1*E1') - (D1*B1' + E1*C1')*A1^(-1) + D1;
    bound(:,:,n) = A*B^(-1)*A' + rrT;
    %     sqrt(bound)
end

%% plot
RMSE1 = RMSE(1,1,:); RMSE1 = RMSE1(:)';
% RMSE2 = RMSE(2,2,:); RMSE2 = RMSE2(:)';
% RMSE3 = RMSE(3,3,:); RMSE3 = RMSE3(:)';
bound1 = bound(1,1,:); bound1 = bound1(:)';
% bound2 = bound(2,2,:); bound2 = bound2(:)';
% bound3 = bound(3,3,:); bound3 = bound3(:)';
close all
figure(1),
hold on
% semilogy(hCandidate,RMSE1,hCandidate,sqrt(bound1))
% semilogy(rouCandidate,RMSE2,rouCandidate,sqrt(bound2))
% semilogy(rouCandidate,RMSE3,rouCandidate,sqrt(bound3))
semilogy(hCandidate,RMSE1,'LineWidth',1.2,'Marker','*')
semilogy(hCandidate,sqrt(bound1),'LineWidth',2)
legend('RMSE of $\hat{\theta}_1 - \psi_1$','Bound of $\theta_1 - \psi_1$',...
     'FontSize',12,'Interpreter','latex')
xlabel('h',FontSize=14)
ylabel('Estimation Error',FontSize=14)



