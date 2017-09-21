%
% myShelf
% 
% Calculate digital shelf filter coefficients
%    [b,a]   = myShelf(N,wc,g,type),  - OR -
%    [z,p,k] = myShelf(N,wc,g,type)
%
% Inputs
%    N    - filter order
%    wc   - normalized digital cutoff frequency between 0 and 1.  If
%           lowpass or highpass this should be a scalar, but if bandpass or
%           bandstop this should have two elements for the two cutoffs
%    g    - DC gain dB
%    type - string specifying the filter type as 'low', 'high', 'bandpass',
%           or 'bandstop'.  Leave this out to assume low pass.
%
% Outputs
%    If 2 outputs requested, [b,a] vectors returned in transfer 
%       function form
%    If 3 outputs requested, [z,p,k] vectors and gain returned in 
%       zero-pole-gain form
%
% References: 
%   "Parametric Higher Order Shelving Filters"
%   Martin Holters and Udo Zolzer
%   14th European Signal Processing Conference (EUSIPCO 2006)
%   Florence, Italy, September 4-8, 2006
%
% Written by Phil Townsend (jptown0@gmail.com) 11/21/2014 

function [b,a,k] = myShelf(N,wc,g,type)
    wc    = tan(pi/2*wc);         % prewarp cutoff frequency
    n     = (1:N)';               % iterator
    alpha = pi/2*(1-(2*n-1)/N);   % complex exponential arg
    wo2   = prod(wc);             % cutoff (squared)
    B     = abs(diff(wc));        % bandwidth   
    g     = 10^(g/20);            % convert gain to linear  

    % Design an analog prototype first, then use the Bilinear Transform to
    % map to digital coefficients
    if (nargin<4) || strcmp(type,'low')
        pa = -wc*exp(j*alpha)*(g^(-1/2/N)); % analog poles
        za = -wc*exp(j*alpha)*(g^(1/2/N));  % analog zeros
        ka = 1;                             % analog gain
    elseif strcmp(type,'high')
        pa = -wc*exp(-j*alpha)*(g^(1/2/N));   % analog poles
        za = -wc*exp(-j*alpha)*(g^(-1/2/N));  % analog zeros
        ka = g;                               % analog gain
    elseif strcmp(type,'bandpass')
        if length(wc) ~= 2
           error('wc must have two elements') 
        end
        pa = myQuadratic(ones(N,1), B*exp(j*alpha), wo2*ones(N,1)); % analog poles
        pa = sortPoles(pa);
        za = myQuadratic(ones(N,1), B*exp(j*alpha)*(g^(1/N)), wo2*ones(N,1)); % analog zeros
        za = sortPoles(za);
        ka = 1; % analog gain
    elseif strcmp(type,'stop')
        if length(wc) ~= 2
           error('wc must have two elements') 
        end        
        pa = myQuadratic(exp(j*alpha), B*ones(N,1), wo2*exp(j*alpha)); % analog zeros
        pa = sortPoles(pa);
        za = myQuadratic(exp(j*alpha)*(g^(1/N)), B*ones(N,1), wo2*exp(j*alpha)*(g^(1/N))); % analog zeros
        za = sortPoles(za);
        ka = g; % analog gain
    else
        error('Unknown type %s', type)
    end        
    pd = (1+pa)./(1-pa);  % pole bilinear transform
    zd = (1+za)./(1-za);  % zero bilinear transform
    zd = [zd; -ones(length(pd)-length(zd),1)];  % Add extra zeros at -1
    kd = real(ka*prod(1-za)./prod(1-pa));       % gain bilinear transform
    if (nargout == 3)     % pass [z,p,k] form to output
        b = zd;  a = pd;  k = kd;
    else  % two outputs, pass [b,a] to output
        b  = real(kd*poly(zd)); % TF num from polynomial mult
        a  = real(poly(pd));    % TF denom from polynomial mult
    end
end

% Quadratic Formula
% Returns x = (-b +/- sqrt(b^2-4ac)) / 2a for vectors a, b, and c
% Inefficient Matlab style but helps ensure poles kept in order
function x = myQuadratic(a,b,c)
    x = zeros(2*length(a),1);
    for k = 1:length(a)
        x(2*k-1) = (-b(k) + sqrt(b(k).^2-4*a(k)*c(k))) / (2*a(k));
        x(2*k)   = (-b(k) - sqrt(b(k).^2-4*a(k)*c(k))) / (2*a(k));
    end
end

% Get bandpass poles in strictly complex conjugate pairs
function y = sortPoles(x)
    y = x;
    for n = 1:4:length(x)-3
        y(n+(0:3)') = x(n+[0;2;1;3]);
    end
end
