%
% cardioidLinearArray.m
% 
% DESCRIPTION
%   This script computes and plots the beampattern of a uniform linear
%   array (ULA) of M directional elements dm meters apart steered to an
%   angle phi from endfire.  The elements supported are modelled as first
%   order differential microphones of two elements spaced dc meters apart
%   and oriented broadside, like this:
%
%    ^       ^       ^       ^
%   / \     / \     / \     / \
%    |-------|-------|-------|
%    |       |       |       |
%
%   The combination of cardioid-type elements into
%   a summing array is a nicely complementary combination since
%   differential mics have good low freuqency directivity but suffer SNR
%   loss while summing arrays have good high frequency directivity and
%   improve SNR.
%
% PARAMETERS
%   f       - frequency of interest in Hz
%   c       - speed of sound in m/s.  Typically 340-343 m/s
%   dm      - linear array spacing in meters
%   dc      - cardioid internal element spacing in meters
%   M       - number of directional microphones in the linear array
%   dbLim   - polar axis lower limit in dB
%   eltType - one of 'omni', 'cardioid', 'hypercard', or 'fig8'
%
% REFERENCES
%   "Superdirectional Microphone Arrays" by Gary Elko, Ch 10 of "Acoustic
%   Signal Processing for Telecommunications", pg 181-237
%   "Microphone Arrays: A Tutorial" by Ian McCowan
%
% Written by Phil Townsend (jptown0@gmail.com) Sep 12, 2017

%% Params
theta   = linspace(0,2*pi,1e3)'; % angle axis (radians)
phi     = pi/3;                  % steering angle (radians; broadside is pi/2)
f       = 1e3;                   % frequency (Hz)
c       = 343;                   % speed of sound (m/s)
dm      = .035;                  % intermic spacing (m)
dc      = .021;                  % cardioid interelement spacing (m)
M       = 8;                     % number of directional mics in the array
m       = dm*(1:M);              % mic positions (m)
dbLim   = 20;                    % dB limit for polar plots
eltType = 'cardioid';            % 'omni', 'cardioid', 'hypercard', 'fig8'

%% Calculate
% Select element response
switch eltType
    case 'omni',      wfun = @(theta)(ones(size(theta)));
    case 'cardioid',  wfun = @(theta)(1-exp(j*pi*f*dc/c*(1+sin(theta))));
    case 'hypercard', wfun = @(theta)(1-exp(j/2*pi*f*dc/c*(1+3*sin(theta))));
    case 'fig8',      wfun = @(theta)(1-exp(j*2*pi*f*dc/c*sin(theta)));
    otherwise,        error('Unknown element type')
end

% Array Response
B = sum(wfun(theta).*exp(-j*2*pi*f/c*m.*(cos(theta)-cos(phi))),2);
B = 20*log10(abs(B));
B = B - max(B) + dbLim;
B(B<0) = 0;
figure, subplot(1,2,1), h = polar(theta,B);  set(h,'linewidth',2)
title(sprintf('ULA Beampattern, f = %g Hz, M = %g, \\phi = %g\\circ',...
    f, length(m), phi*180/pi))

% Element Response
Bw = 20*log10(abs(wfun(theta)));
Bw = Bw - max(Bw) + dbLim;
Bw(Bw<0) = 0;
subplot(1,2,2), h = polar(theta, Bw);  set(h,'linewidth',2)
title('Element Response')

