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
    [x, order] = sort(x);
    y = y(order);

    % Collapse repeated abscissas by averaging their ordinates so every
    % interpolation interval has a nonzero width.
    [xUnique, ~, group] = unique(x);
    if numel(xUnique) < 2
        error('linearinterp:InsufficientUniqueData', ...
            'x must contain at least two distinct values.');
    end
    if numel(xUnique) < numel(x)
        y = accumarray(group, y, [], @mean);
        x = xUnique;
    end

    yq = interp1(x, y, xq(:), 'linear', 'extrap');
    yq = reshape(yq, outputSize);
end