% LINEARINTERP Linearly interpolate and extrapolate data.
%
% -- Usage --
% yq = LINEARINTERP(x,y,xq)
%
% -- Input --
% x,y = vectors describing a function y = f(x);
%       x does not need to be sorted.
% xq  = query points at which to evalulate the output
%
% -- Output --
% yq  = vector of function values at the query points: yq = f(xq) 
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

function yq = linearinterp(x, y, xq)
  
    % Initialize
    x = x(:);    % force column vectors
    y = y(:);
    
    size_original_xq = size(xq);
    xq = xq(:);
    yq = zeros(size_original_xq);
    
    % Sort x to ensure it is monotonic for interpolation.
    % (noise may corrupt underlying monotonic function).
    % Strictly speaking, this alters the function y = f(x),
    % but the alterations are within the noise.
    [x, IX] = sort(x);
    y       = y(IX);
    
    % Indicies to three different types of x query points.
    lo  = xq <= x(1);    % logical indicies to points off the low end
    hi  = xq >= x(end);  % logical indicies to points off the high end
    mid = find(x(1) < xq & xq < x(end));  % indicies to mid points
    
    % 1. xq < x(1): extrapolate off low end of table
    if any(lo)
        dx     = x(2) - x(1);
        dy     = y(2) - y(1);
        yq(lo) = (xq(lo)-x(1))*dy/dx + y(1);
    end
    
    % 2. xq > x(end): extrapolate off high end of table
    if any(hi)
        dx     = x(end) - x(end-1);
        dy     = y(end) - y(end-1);
        yq(hi) = (xq(hi)-x(end))*dy/dx + y(end);
    end
    
    % 3. x(1) < xq < x(end): linearly interpolate
    if any(mid)
        % Find interval index for each xq
        % bin(i) = k such that x(k) <= xq(i) < x(k+1)
        [~, ~, loc] = histcounts(xq(mid), x);
    
        % Fix edge case loc = length(x)
        loc(loc == 0) = 1;
        loc(loc == numel(x)) = numel(x)-1;
    
        x1 = x(loc);
        x2 = x(loc+1);
        y1 = y(loc);
        y2 = y(loc+1);
    
        m       = (y2 - y1) ./ (x2 - x1);
        yq(mid) = y1 + m .* (xq(mid) - x1);
    end