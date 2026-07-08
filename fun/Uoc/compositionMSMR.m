function [Z, dZ, xj, Uj, dZj] = compositionMSMR(Uvect, params, TdegC)
  
% Constants.
R = 8.3144598;  % Molar gas constant [J/mol K]
F = 96485.3329; % Faraday constant [C/mol]

U0 = params.U0(:);
X = params.X(:);
w = params.w(:);

f = F/R/(TdegC + 273.15);
Uvect = Uvect(:).';

gj = exp(f*(Uvect - U0)./w);
xj = X./(1 + gj);
Uj = U0 + w./f.*log((X - xj)./xj);

Z = sum(xj, 1);

dZj = (X ./ w) .* gj ./ (1 + gj) .^ 2;
dZj(isnan(dZj)) = 0; % !!! when gj->Inf, (gj./(1+gj).^2)->0
dZ = f * sum(dZj, 1);

Z = Z(:);
dZ = dZ(:);

