% COMPOSITIONMSMR Calculate MSMR composition and differential composition.
%
% -- Usage --
% [Z,dZ,xj,Uj,dZj] = compositionMSMR(Uvect,params,TdegC)
%
% -- Description --
% Calculates the composition of an electrode using a multi-site
% multi-reaction (MSMR) model.
%
% Each MSMR gallery j is defined by:
%   U0(j) = equilibrium potential of gallery j [V]
%   X(j)  = capacity fraction of gallery j [-]
%   w(j)  = transition-width parameter
%
% For a given potential U, the function calculates the occupation of each
% gallery xj, the corresponding gallery potential Uj, and the differential
% contribution dZj/dU. The total electrode composition is obtained as the
% sum of the individual gallery occupations.
%
% -- Inputs --
% Uvect  = Electrode potential vector [V].
% params = Structure containing the MSMR parameters:
%          .U0       Gallery equilibrium potentials [V].
%          .X        Gallery capacity fractions [-].
%          .w        Gallery transition-width parameters.
%
% TdegC  = Temperature [degC].
%
% -- Outputs --
% Z      = Total electrode composition as a function of potential [-].
% dZ     = Positive differential magnitude -dZ/dU [1/V], N-by-1.
% xj     = Occupation of each MSMR gallery, J-by-N, dimensionless.
% Uj     = Supplied potentials repeated for each gallery, J-by-N [V].
% dZj    = Positive gallery derivative magnitudes divided by f=F/(R*T),
%          J-by-N. Multiply by f to obtain |dx_j/dU| [1/V].
%
% -- Model --
% The total composition is calculated as:
%   Z(U) = sum_j x_j(U)
% where x_j is the contribution of gallery j.
% The signed differential composition is:
%   dZ/dU = sum_j dx_j/dU
% The returned dZ is its negative, i.e. the positive magnitude.
% The calculation assumes positive gallery fractions X and positive
% transition-width parameters w.

function [Z, dZ, xj, Uj, dZj] = compositionMSMR(Uvect, params, TdegC)

validateattributes(Uvect, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
validateattributes(TdegC, {'numeric'}, {'scalar', 'real', 'finite', '>', -273.15});
% Physical constants
R = 8.3144598;       % Molar gas constant [J/(mol K)]
F = 96485.3329;      % Faraday constant [C/mol]

% MSMR parameters
U0 = params.U0(:);
X  = params.X(:);
w  = params.w(:);
validateattributes(U0, {'numeric'}, {'vector', 'real', 'finite', 'nonempty'});
validateattributes(X, {'numeric'}, {'vector', 'real', 'finite', 'positive'});
validateattributes(w, {'numeric'}, {'vector', 'real', 'finite', 'positive'});
if numel(U0) ~= numel(X) || numel(U0) ~= numel(w)
    error('compositionMSMR:ParameterSize', 'U0, X, and w must have matching lengths.');
end

% Temperature-dependent factor
f = F/R/(TdegC + 273.15);

% Ensure potential is a row vector
Uvect = Uvect(:).';

% A sign-aware logistic avoids overflow and cancellation in saturated sites.
argument = f*(Uvect - U0)./w;
tail = exp(-abs(argument));
occupation = tail./(1 + tail);
negative = argument < 0;
occupation(negative) = 1./(1 + tail(negative));
xj = X.*occupation;
% These are exactly the input potentials, even where occupation rounds to 0/1.
Uj = repmat(Uvect, numel(U0), 1);

% Total composition
Z = sum(xj,1);

% Differential contribution of each gallery
dZj = (X./w) .* tail ./ (1 + tail).^2;

% Total differential composition
dZ = f*sum(dZj,1);

% Return column vectors
Z  = Z(:);
dZ = dZ(:);

