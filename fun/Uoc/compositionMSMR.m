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
% dZ     = Differential electrode composition dZ/dU [1/V].
% xj     = Occupation of each MSMR gallery as a function of potential.
% Uj     = Potential corresponding to the occupation of each gallery [V].
% dZj    = Differential contribution of each MSMR gallery to dZ/dU [1/V].
%
% -- Model --
% The total composition is calculated as:
%   Z(U) = sum_j x_j(U)
% where x_j is the contribution of gallery j.
% The differential composition is:
%   dZ/dU = sum_j dx_j/dU
% The calculation assumes positive gallery fractions X and positive
% transition-width parameters w.

function [Z, dZ, xj, Uj, dZj] = compositionMSMR(Uvect, params, TdegC)
  
% Physical constants
R = 8.3144598;       % Molar gas constant [J/(mol K)]
F = 96485.3329;      % Faraday constant [C/mol]

% MSMR parameters
U0 = params.U0(:);
X  = params.X(:);
w  = params.w(:);

% Temperature-dependent factor
f = F/R/(TdegC + 273.15);

% Ensure potential is a row vector
Uvect = Uvect(:).';

% Gallery occupation
gj = exp(f*(Uvect - U0)./w);
xj = X./(1 + gj);

% Potential corresponding to each gallery occupation
Uj = U0 + w./f .* log((X - xj)./xj);

% Total composition
Z = sum(xj,1);

% Differential contribution of each gallery
dZj = (X./w) .* gj ./ (1 + gj).^2;

% Numerical limit when gj -> Inf
dZj(isnan(dZj)) = 0;

% Total differential composition
dZ = f*sum(dZj,1);

% Return column vectors
Z  = Z(:);
dZ = dZ(:);

