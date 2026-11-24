% LINEARINTERP Linearly interpolate and extrapolate data.
%
% -- Usage --
% yq = linearinterp(x,y,xq)
%
% -- Input --
% x,y = vectors describing a function y = f(x);
%       x does not need to be sorted.
% xq  = query points at which to evaluate the output; any array shape.
%
% -- Output --
% yq  = function values at the query points, with the same shape as xq.
% Duplicate x values are combined by averaging their y values. At least
% two distinct x values are required; extrapolation is linear.
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

% Modified by Sergio Fernandez Gonzalez: input validation, duplicate-sample
% averaging, query-shape preservation, and documentation updates.
%
function yq = linearinterp(x, y, xq)

    validateattributes(x, {'numeric'}, {'vector','real','finite','nonempty'}, mfilename, 'x');
    validateattributes(y, {'numeric'}, {'vector','real','finite','nonempty'}, mfilename, 'y');
    validateattributes(xq, {'numeric'}, {'real','finite'}, mfilename, 'xq');
    if numel(x) ~= numel(y)
        error('linearinterp:InputSizeMismatch', 'x and y must contain the same number of elements.');
    end
    if numel(x) < 2
        error('linearinterp:InsufficientData', 'x and y must contain at least two samples.');
    end

    outputSize = size(xq);
    x = x(:);
    y = y(:);
    % Already increasing data need no sorting or duplicate grouping.
    if any(x(2:end) < x(1:end-1))
        [x, order] = sort(x);
        y = y(order);
    end

    % After sorting, repeated abscissas are adjacent. Average their ordinates
    % only when needed so every interpolation interval has a nonzero width.
    if any(x(2:end) == x(1:end-1))
        [xUnique, ~, group] = unique(x);
        if numel(xUnique) < 2
            error('linearinterp:InsufficientUniqueData', ...
                'x must contain at least two distinct values.');
        end
        y = accumarray(group, y, [], @mean);
        x = xUnique;
    end

    yq = interp1(x, y, xq(:), 'linear', 'extrap');
    yq = reshape(yq, outputSize);
end
