%% Setup Parameters
clear; clc

%% Setup Experiment

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
param_th.h = 1;
param_th.Sigma = 0.1*eye(dimension);

% Storing Results
ps = zeros(dimension,n_trial);
% x = zeros(3,param_ps.xSample,n_trial);
est_th = zeros(dimension,n_trial);
theta0 = zeros(dimension,n_trial);

% Simulation Variable Setup
NOSCandidate = (1:1:20);
nCandidate = length(NOSCandidate);
RMSE = zeros(dimension,dimension, nCandidate);
bound = zeros(dimension,dimension, nCandidate);
% CRB = zeros(3,3, nCandidate);

for n = 1:nCandidate
    n
    nSamples = NOSCandidate(n);
    param_ps.xSample = nSamples;
    param_th.xSample = nSamples;
    for m = 1:n_trial
        %% Data Generation
        % Data of Prior
        ps(:,m) = mvnrnd(param_ps.mu,param_ps.sigma);
        % Data of x distribution
        x = mvnrnd(param_ps.h*ps(:,m),param_ps.Sigma,nSamples)';
        % Pseudotrue
        theta0(:,m) = (nSamples*param_th.h'*param_th.Sigma^(-1)*param_th.h + param_th.sigma^(-1))^(-1)...
            *(nSamples*param_th.h'*param_th.Sigma^(-1)*param_ps.h*ps(:,m)+ param_th.sigma^(-1)*param_th.mu);
        % Estimation
        est_th(:,m) = (nSamples*param_th.h'*param_th.Sigma^(-1)*param_th.h + param_th.sigma^(-1))^(-1)...
            *(param_th.h'*param_th.Sigma^(-1)*sum(x,2)+ param_th.sigma^(-1)*param_th.mu); % Implement MAP
        %         mean_est_th(:,m) = mean(est_th(1,1:m)); % mean of estimate convergence
    end

    %% Error Record
    error = est_th - theta0;
    covar = zeros(dimension,dimension);
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
    bound(:,:,n) = A*B^(-1)*A' ;%+ rrT;
    b(:,:,n) = B^(-1);
    %     sqrt(bound)
end

%% plot
RMSE1 = RMSE(1,1,:); RMSE1 = RMSE1(:)';
RMSE2 = RMSE(2,2,:); RMSE2 = RMSE2(:)';
RMSE3 = RMSE(3,3,:); RMSE3 = RMSE3(:)';
bound1 = bound(1,1,:); bound1 = bound1(:)';
bound2 = bound(2,2,:); bound2 = bound2(:)';
bound3 = bound(3,3,:); bound3 = bound3(:)';
b1 = b(1,1,:); b1 = b1(:)';

close all
figure(1),
hold on
semilogy(NOSCandidate,RMSE1,'LineWidth',1.2,'Marker','*')
semilogy(NOSCandidate,sqrt(bound1),'LineWidth',2)
% semilogy(NOSCandidate,RMSE2,NOSCandidate,sqrt(bound2))
% semilogy(NOSCandidate,RMSE3,NOSCandidate,sqrt(bound3))
legend('RMSE of $\hat{\theta}_1 - \psi_1$','Bound of $\theta_1 - \psi_1$',...
     'FontSize',12,'Interpreter','latex')
%     'RMSE of $\hat{\theta}_2 - \psi_2$','MBCRB of $\theta_2 - \psi_1$',...
%     'RMSE of $\hat{\theta}_3 - \psi_3$','MBCRB of $\theta_3$ - \psi_1',
       
xlabel('Number of Samples',FontSize=14)
ylabel('Estimation Error',FontSize=14)

close all
figure(2)
hold on
semilogy(NOSCandidate,RMSE1,'LineWidth',1.2,'Marker','*')
semilogy(NOSCandidate,sqrt(bound1),'LineWidth',2)
semilogy(NOSCandidate,sqrt(b1),'LineWidth',2)
legend('RMSE of $\hat{\theta}_1$','MBCRB of $\theta_1$',...
    'BCRB of $\hat{\theta}_1$','FontSize',12,'Interpreter','latex')
xlabel('Number of Samples',FontSize=14)
ylabel('Estimation Error',FontSize=14)



