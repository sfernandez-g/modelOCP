% SMOOTHDIFF Smooth a voltage-v-soc curve and estimate the differential
% capacity.
%
% -- Usage --
% [zout,vout,dz] = smoothdiff(z,v)
% [zout,vout,dz] = smoothdiff(z,v,dv)
% [zout,vout,dz] = smoothdiff(z,v,dv,tol)
% [xout,vout,dx] = smoothdiff(x,v,dv,tol,n)
%
% -- Input --
% z    = soc vector [unitless]
% v    = voltage vector [V]
% dv   = voltage resolution (default 1e-3) [V]
% tol  = edge tolerance (default 5e-5) [V]
% n    = number of points used to resample the x-v curve
%        (default = number of input points) [unitless]
%
% -- Output --
% zout = smoothed soc vector [unitless]
% vout = smoothed voltage vector [V]
% dz   = differential capacity vector, dz/dv [1/V]
%
% Copyright (©) 2024 The Regents of the University of Colorado, a body
% corporate. Created by Gregory L. Plett and M. Scott Trimboli of the
% University of Colorado Colorado Springs (UCCS). This work is licensed
% under a Creative Commons "Attribution-ShareAlike 4.0 International" Intl.
% License. https://creativecommons.org/licenses/by-sa/4.0/
% This code is provided as a supplement to: Gregory L. Plett and M. Scott
% Trimboli, "Battery Management Systems, Volume III, Physics-Based
% Methods," Artech House, 2024. It is provided "as is", without express or
% implied warranty. Attribution should be given by citing: Gregory L. Plett
% and M. Scott Trimboli, Battery Management Systems, Volume III:
% Physics-Based Methods, Artech House, 2024.
%
% Modified by Fernandez Gonzalez Sergio:
% Added "x array monotonical" and some comments

function [xout, vout, dx] = smoothdiff(x, v, dv, tol, n)

validateattributes(x, {'numeric'}, {'vector','real','finite'}, mfilename, 'x');
validateattributes(v, {'numeric'}, {'vector','real','finite','numel',numel(x)}, mfilename, 'v');
if numel(x) < 2
    error('smoothdiff:InsufficientData', 'At least two paired input samples are required.');
end
x = x(:);
v = v(:);

% Remove consecutive duplicate x values.
vec = [true; diff(x) ~= 0];
x = x(vec);
v = v(vec);

% Initialization.
if nargin < 3 || isempty(dv)
    dv = 1e-3;
end
if nargin < 4 || isempty(tol)
    tol = 5e-5;
end
if nargin < 5 || isempty(n)
    n = numel(v);
end
validateattributes(dv, {'numeric'}, {'scalar','positive','finite'}, mfilename, 'dv');
validateattributes(tol, {'numeric'}, {'scalar','nonnegative','finite'}, mfilename, 'tol');
validateattributes(n, {'numeric'}, {'scalar','integer','finite','>=',2}, mfilename, 'n');

% Define the new strictly monotonic x grid
xmin = min(x);
xmax = max(x);
xgrid = linspace(xmin, xmax, n).';
v = linearinterp(x, v, xgrid);
x = xgrid;

% Define voltage points.
voltageMin = min(v);
voltageMax = max(v);
vmin = voltageMin - tol;
vmax = voltageMax + tol;
vout = (vmin:dv:vmax)';

% Approximate differential capacity dz/dv.
Nv = histcounts(v, numel(vout));
dx = (xmax-xmin)/n*Nv/dv;

% Interpolate output soc.
xout = linearinterp(v, x, vout);

% Force the boundaries to agree with the original data-set.
[~, ind0] = min(xout);
xout(ind0) = xmin;
[~, ind1] = max(xout);
xout(ind1) = xmax;
[~, ind0] = min(vout);
vout(ind0) = voltageMin;
[~, ind1] = max(vout);
vout(ind1) = voltageMax;

% Force column vector
xout = xout(:);
vout = vout(:);

kernel = ones(2, 1) / 2;
dx = conv(dx, kernel, 'same');
dx = dx(:);

end
