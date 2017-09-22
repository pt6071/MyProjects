%
% scatterPlot.m
%
% Plots the pressure field for a plane wave propagating up the y axis
% scattered by either a sphere or an infinitely-tall rigid cylinder. This
% is a great way of visualizing the interferece pattern from the
% reflections as well as the acoustic shadows behind the obstacle and
% comparing how the sphere and cylinder each affect the field.
%
% Written by Phil Townsend (jptown0@gmail.com) 4/25/17

%% Parameters
R  = .01;       % radius of the cylinder
c  = 343;       % speed of sound in m/s
f  = 20e3;      % frequency of interest
k  = 2*pi*f/c;  % wavenumber
N  = 10;        % inf sum truncation
ax = linspace(-.05,.05,500)';  % linear axis for plots

% discfn plots a disc representing the cylinder on pressure field plots
discfn = @(R)(patch([linspace(-R,R,100) -linspace(-R,R,100)], ...
    [-sqrt(R^2-linspace(-R,R,100).^2) sqrt(R^2-linspace(-R,R,100).^2)], [1 1 1]));

%% Solve for Pressure Fields
P = cylindricalSolver(N,k,R,ax);  % faster as a function
figure, hold on, imagesc(ax,ax,20*log10(abs(P)).'), axis tight square
hd = colorbar;  set(get(hd,'ylabel'),'String','Power (dB)');  discfn(R)
title(sprintf('Cylindrical Scattering Pressure Field Magnitude\nf = %g Hz, R = %g m', f, R))
xlabel('x (m)'), ylabel('y (m)'), colormap(jet)

P = sphericalSolver(N, ax, k, R);  % faster as a function
figure, imagesc(ax,ax,20*log10(abs(P))), axis square
title(sprintf('Spherical Scattering Pressure Field Magnitude\nf = %g Hz, R = %g m', f, R))
xlabel('x (m)'), ylabel('y (m)'), colormap(jet), set(gca,'YDir','normal')
hd = colorbar;  set(get(hd,'ylabel'),'String','Power (dB)');  discfn(R)
