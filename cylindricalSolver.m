%
% cylindricalSolver.m
%
% Description
%    Solver function for computing the 2D pressure field at wavenumber k
%    scattered by a rigid cylinder of radius R and infinite height on the
%    square domain defined at the points in ax truncating the infinite
%    summation to N terms.
%
% References
%   Section 1.2.3 of "Analytical Techniques for Acoustic Scattering by
%   Arrays of Cylinders" by Nikolaos Tymis, Doctoral Thesis
%
% Written by Phil Townsend (jptown0@gmail.com) 5/15/17

function P = cylindricalSolver(N,k,R,ax)
    [xx,yy]  = meshgrid(ax);       % mesh coords for 2D plot
    rr       = sqrt(xx.^2+yy.^2);  % radial distance for calculations
    rr(rr<R) = nan;                % avoid calculations inside cylinder
    tt       = atan2(yy,xx);       % angles
    n        = (-N:N)';            % n axis
    
    % cat is used to take advantage of symmetry
    tmp  = exp(pi/2*j*n).*(besselj(n-1,k*R)-besselj(n+1,k*R)) ./ ...
        (besselh(n-1,1,k*R)-besselh(n+1,1,k*R));
    tmp  = permute(tmp(:,ones(length(ax),1),ones(length(ax),1)),[3 2 1]);
    n    = (0:N)';   % setup for rest of calc's
    nn   = permute(n(:,ones(length(ax),1),ones(length(ax),1)),[3 2 1]);
    ejnt = exp(j*nn.*tt(:,:,ones(size(n))));  % calculate first half
    ejnt = cat(3,1./ejnt(:,:,end:-1:2),ejnt); % use reciprocal symmetry
    Hnkr = besselh(nn,1,k*rr(:,:,ones(size(n)))); % calculate first half
    Hnkr = cat(3,(-1).^nn(:,:,end:-1:2).*Hnkr(:,:,end:-1:2),Hnkr); % symmetry
    P    = exp(j*k*xx) - sum(tmp.*Hnkr.*ejnt,3);  % incident and scattered waves
end
