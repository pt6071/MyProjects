%
% blms.m
%
% Implements a block (frequency domain) LMS filter.
%
% Inputs
%    u - Input signal
%    d - Desired signal
%    M - order of the filter (power of 2 is good for FFT's)
%    alpha - adaptation constant (internally scaled by signal power)
%    gamma - forgetting factor of power estimates (for nonstationary sigs)
%
% Outputs
%    y - Filtered input signal
%    W - frequency domain filter taps every iteration
%
% Based on Sec 10.2 of "Adaptive Filter Theory, 3rd Ed" by Haykin
% Written by Phil Townsend (jptown0@gmail.com) 10/25/2012

function [y, e, w] = blms(u, d, M, alpha, gamma)
    % Initialization
    nBlocks = floor(length(d)/M)-1; % total number of blocks to iterate
    y = zeros(size(d)); % overall output vector
    e = zeros(size(d)); % error signal
    w = zeros(2*M,1);   % tap weights each block
    P = zeros(2*M,1);   % power estimate each block (for scaling)
    uWin = u(1:M);      % store first block of input data
    
    % Block by block filtering and adaptation
    for k = 1:nBlocks     % consider the 0th block an init step
        uWinPrev = uWin;  % hand off previous window of input
        uWin = u(k*M+1:(k+1)*M);  % store current window of input
        U    = diag(fft([uWinPrev; uWin]));  % input in freq domain
        yWin = ifft(U*w); % filter input and take into time domain
        y(k*M+1:(k+1)*M) = yWin(end-M+1:end); % keep last M elts of yWin
        dWin = d(k*M+1:(k+1)*M);  % store current window of reference
        eWin = dWin - yWin(end-M+1:end);  % window of error signal 
        e(k*M+1:(k+1)*M) = eWin;  % retain error signal for analysis
        E    = fft([zeros(M,1); eWin]);  % error sig in freq domain
        P    = (1-gamma)*P + gamma*abs(diag(U)).^2; % power estimate update
        D    = diag(1./P);    % make diag matrix for correlation
        phi  = ifft(D*U'*E);  % correlation (scaled by power)
        phi  = phi(1:M);      % only want first M samples of this
        w    = w + alpha*fft([phi; zeros(M,1)]); % tap update
    end    
end

