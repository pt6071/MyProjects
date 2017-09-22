%
% nlms.m
%
% Basic implementation of a leaky normalized LMS adaptive filter
% 
% Inputs
%    x    - Input signal
%    d    - Desired signal
%    mu   - adaptive step size (0 < mu < 2)
%    ord  - Order of the filter
%    beta - Forgetting factor (0 < beta < 1)
%
% Outputs
%    y    - filtered output signal
%    e    - error signal for all input samples
%    w    - final tap weight vector
%
% Written by Phil Townsend (jptown0@gmail.com) 10/8/2012

function [y, e, w] = nlms(x, d, mu, ord, beta)
    % Initializations
    y = zeros(length(x),1);   % filter output
    e = zeros(length(x),1);   % error signal
    w = zeros(1,ord);         % initialize taps
    
    % Filter and updates
    for n = ord:length(x)
        xWin = x(n-ord+1:n);
        y(n) = w*xWin;
        e(n) = d(n)-y(n);
        w = beta*w + mu*xWin'*e(n)/(xWin'*xWin);
    end
end