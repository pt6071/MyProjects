%
% blmsdemo.m
%
% A quick demonstration of system identification comparing the Normalized
% LMS and Block LMS algorithms.  NLMS filters are very simple, robust, and
% reliable while BLMS filters are more complicated but less computationally
% expensive for long filter lengths thanks to internal use of FFT's.  Here
% we adapt a fairly long filter to show that BLMS can yield a lower steady
% state mean squared error (MSE) even if NLMS can initially learn more
% quickly by plotting the ensemble-averaged learning curves for both.
%
% Written by Phil Townsend (jptown0@gmail.com) 9/21/17

%% Params
Nruns = 500;             % number of runs for ensemble average
alpha = .6;              % BLMS step size
gamma = .5;              % BLMS forgetting factor 
mu    = 1;               % NLMS step size (0 < mu < 2)
beta  = 1;               % NLMS forgetting factor (0 < beta < 1)
M     = 64;              % number of filter taps to estimate
N     = 1024;            % number of samples to process
[b,a] = cheby1(5,10,.5); % system to model
eBlms = zeros(N,Nruns);  % BLMS learning curves for ensemble average 
eNlms = zeros(N,Nruns);  % NLMS learning curves for ensemble average

%% Run and Plot
for nn = 1:Nruns
    x = randn(N,1);      % create some Gaussian white noise
    d = filter(b,a,x);   % pass the noise through the system
    [~, eBlms(:,nn), ~] = blms(x, d, M, alpha, gamma);  % apply 
    [~, eNlms(:,nn), ~] = nlms(x, d, mu, M, beta);      % apply
end
figure, plot(1:N,10*log10(mean(eBlms.^2,2)),1:N,10*log10(mean(eNlms.^2,2)))
grid on, axis tight
title('Adaptive Filter Learning Curves')
xlabel('Samples'), ylabel('MSE (dB)')
legend('Block LMS','Normalized LMS')
