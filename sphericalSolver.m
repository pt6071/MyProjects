%
% sphericalSolver.m
%
% Description
%    Solver function for computing the 2D pressure field at wavenumber k
%    scattered by a rigid sphere of radius R and infinite height on the
%    square domain defined at the points in ax truncating the infinite
%    summation to N terms.
%
% References
%    Section 12.3 of "Acoustics: Sound Fields and Transducers" by Beranek
%    and Mellow
%
% Written by Phil Townsend (jptown0@gmail.com) 12/30/2016

function P = sphericalSolver(N, ax, k, R)
    % Setup
    n          = 0:N;               % summation axis
    [xm,ym]    = meshgrid(ax);      % 2D mesh for calculations
    r          = sqrt(xm.^2+ym.^2); % convert to polar (field only)
    cos_t      = ym./r;             % easy evaluation of cosine
    P          = exp(-j*k*ym);       % incoming plane wave
    [xm,ym,nm] = ndgrid(ax,ax,n);   % xy mesh for each n
    rm         = sqrt(xm.^2+ym.^2); % r  mesh for each n
    rm(rm<R)   = nan;               % don't compute inside
    
    % Scattered Field
    jpn  = n.*besselj(n-.5,k*R)-(n+1).*besselj(n+1.5,k*R);
    hpn2 = n.*(besselj(n-.5,k*R)- j*bessely(n-.5,k*R)) ...
        - (n+1).*(besselj(n+1.5,k*R)- j*bessely(n+1.5,k*R));
    jpn  = ((-j).^n).*(2*n+1).*jpn./hpn2;
    hn2  = sqrt(pi/2/k./rm).*(besselj(nm+.5,k*rm)-j*bessely(nm+.5,k*rm));
    Pn = ones(length(ax),length(ax),N+1);
    for nn = 1:N  % can't avoid iteration to get legendre values
        Pn(:,:,nn+1) = myLegendre(nn,cos_t);
    end
    
    % Total Field is plane wave minus scattered field
    P = P - sum(reshape(jpn,1,1,N+1).*hn2.*Pn,3);
end
