% LITHIATIONBALANCEEQUATION Evaluate the MSMR composition residual.
%
% -- Usage --
% y = lithiationBalanceEquation(U,z,params,TdegC)
%
% -- Description --
% Calculates the residual between a target electrode composition z and the
% composition predicted by the MSMR model at a given potential U.
% The function is used together with a root-finding algorithm such as
% fzero to determine the potential U corresponding to a specified
% composition z.
% The equation solved is:
%   Z_MSMR(U) - z = 0
% where Z_MSMR(U) is the electrode composition calculated by
% compositionMSMR.
%
% -- Inputs --
% U      = Electrode potential at which the MSMR composition is evaluated [V].
% z      = Target electrode composition [-].
% params = Structure containing the MSMR model parameters.
% TdegC  = Temperature [degC].
%
% -- Output --
% y      = Composition residual:
%          y = z - Z_MSMR(U)
% A zero of y corresponds to the potential for which the MSMR model reaches
% the target composition z.

function y = lithiationBalanceEquation(U, z, params, TdegC)
Z = compositionMSMR(U, params, TdegC);
y = z - Z;


