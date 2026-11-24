% ROUGHDIFF Approximate differential capacity using first difference.
%
% -- Usage --
% dZ = roughdiff(Z,U)
%
% Z and U are paired vectors with at least two samples. Z is dimensionless
% composition or SOC; U is potential [V] with nonzero consecutive differences.
% dZ is an N-by-1 signed derivative dZ/dU [1/V], averaged from adjacent
% interval slopes. The endpoint slopes are used at the boundaries.
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

% Modified by Sergio Fernandez Gonzalez: documented usage and derivative units.
%
function dZ = roughdiff(Z,U)
  dZ = diff(Z(:))./diff(U(:));
  dZ = ([dZ(1);dZ]+[dZ;dZ(end)])/2; % Avg of fwd/bkwd diffs
end
