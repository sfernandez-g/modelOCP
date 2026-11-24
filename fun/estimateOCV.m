%ESTIMATEOCV Estimate OCV versus SOC from charge and discharge curves.
%
% -- Usage --
% ocvrel = estimateOCV(elec,config,OCPResolution,OCPEstimators)
%
% -- Input --
% elec   = structure with .name, .cal, and .raw. Each .cal entry contains
%          disZ, disV, chgZ, chgV, and TdegC, as returned by processCU.
%          The resistanceBlend method additionally requires .raw entries
%          with script1 and script3 step/voltage data at current transitions.
% config = electrode configuration structure (see fields below).
% OCPResolution = number of points in each interpolation grid.
% OCPEstimators = cell array selecting 'resistanceBlend', 'voltageAverage',
%                 'socAverage', and/or 'diagonalAverage'.
%
% `config` is a struct with the following fields:
%
%   .vmin     = minimum operating voltage [V]
%   .vmax     = maximum operating voltage [V]
%   .datype   = type of diagonal averaging ('useChg','useDis', or 'useAvg')
%   .daxcorrfiltz = function used to limit SOC interval over which xcorr
%                    of charge/discharge dZ-vs-Z curves is evaluated
%   .daxcorrfiltv = function used to limit voltage interval over which xcorr
%                    of charge/discharge dZ-vs-V curves is evaluated
%
%   The .daxcorrfiltz and .daxcorrfiltv functions accept vectors
%   of the electrode SOC or voltage, respectively, and return a vector of
%   logical indices describing where the cross-correlation should be
%   computed. To use multiple distinct intervals and average the xcorr over
%   the intervals, return a cell array with an entry for each interval.
%   These three fields are required only for the diagonalAverage method.
%
% -- Output --
% ocvrel = relative OCV-versus-composition estimates at each temperature
%          (see fields below)
%
% The struct array `ocvrel` has the following fields:
%
%   .elecname  : electrode name supplied in elec.name
%   .TdegC     : temperature [degC]
%   .estimates : struct array with the following fields:
%     .Z       : relative SOC/composition column vector [-]
%     .U       : potential column vector [V]
%     .dZ      : signed derivative dZ/dU [1/V], column vector
%     .name    : name of estimator
%     .elecname: electrode name (same as above)
%     .TdegC   : temperature (same as above) [degC]
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

% Modified by Sergio Fernandez Gonzalez: repository-local interpolation and
% differentiation helpers, diagonal-average processing, and documentation.
%
function ocvrel = estimateOCV(elec, config, OCPResolution, OCPEstimators)

  % Iterate temperatures.
  clear ocvrel;
  for kt = length(elec.cal):-1:1 % Fill backward to allocate the final array size.
    raw = elec.raw(kt);
    cal = elec.cal(kt);
    data = struct;
    data.name  = elec.name;
    data.TdegC = cal.TdegC;

    % Collect unmodified (i.e., "original") dis/charge curves.
    data.orig.disZ = cal.disZ;
    data.orig.disV = cal.disV;
    data.orig.chgZ = cal.chgZ;
    data.orig.chgV = cal.chgV;
    
    % Define standard soc vector and interpolate voltage.
    data.interp.stdZ = linspace(0, 1, OCPResolution).';
    data.interp.disV = linearinterp(data.orig.disZ, data.orig.disV, data.interp.stdZ);
    data.interp.chgV = linearinterp(data.orig.chgZ, data.orig.chgV, data.interp.stdZ);
    
    % Define standard voltage vector and interpolate soc.
    data.interp.stdV = linspace(config.vmin, config.vmax, OCPResolution).';
    data.interp.disZ = linearinterp(data.orig.disV, data.orig.disZ, data.interp.stdV);
    data.interp.chgZ = linearinterp(data.orig.chgV, data.orig.chgZ, data.interp.stdV);
    
    % Store raw script data for use by resistanceBlend estimator.
    data.raw = raw;
    
    % Generate OCP estimates from each estimator.
    clear estimates;
    for ke = length(OCPEstimators):-1:1
      estName = OCPEstimators{ke};

      % Generate OCP estimate.
      estFnName = ['OCP_' estName];
      [Z,U] = feval(estFnName, data, config);

      % Restrict the estimate to the configured operating voltage grid.
      Z = linearinterp(U, Z, data.interp.stdV);
      U = data.interp.stdV;

      % Calculate differential capacity
      dZ = roughdiff(Z,U);

      estimates(ke) = struct('Z',Z,'U',U,'dZ',dZ, ...
                             'name',estName, ...
                             'elecname',data.name, ...
                             'TdegC',data.TdegC);  % save estimate
    end

    % Store estimates at this temperature.
    ocprel_kt = struct;
    ocprel_kt.elecname = data.name;
    ocprel_kt.TdegC = data.TdegC;
    ocprel_kt.estimates = estimates;
    ocvrel(kt) = ocprel_kt;
    
  end  
end
  
  
% OCP Estimators (four methods) -------------------------------------------

function [Z, U] = OCP_resistanceBlend(data, ~)
%OCP_RESISTANCEBLEND Blend branches after compensating current-step voltage jumps.
%   [Z, U] = OCP_resistanceBlend(data, config) uses data.raw.script1/script3,
%   data.orig branch vectors, and data.interp.stdZ. config is unused.
%   Z is the dimensionless SOC column grid; U is the estimated voltage [V].
  
  rawdata = data.raw;
  disZ = data.orig.disZ;
  disV = data.orig.disV;
  chgZ = data.orig.chgZ;
  chgV = data.orig.chgV;
  stdZ = data.interp.stdZ;
  
  % Compute the instant voltage jumps
  indD  = find(rawdata.script1.step == 2); % discharge data (script 1)
  indC  = find(rawdata.script3.step == 2); % charge data (script 3)
  IR1Da = rawdata.script1.voltage(indD(1)-1) - ...
        rawdata.script1.voltage(indD(1));
  IR2Da = rawdata.script1.voltage(indD(end)+1) - ...
        rawdata.script1.voltage(indD(end));
  IR1Ca = rawdata.script3.voltage(indC(1)) - ...
        rawdata.script3.voltage(indC(1)-1);
  IR2Ca = rawdata.script3.voltage(indC(end)) - ...
        rawdata.script3.voltage(indC(end)+1);
  IR1D = min(IR1Da,2*IR2Ca); IR2D = min(IR2Da,2*IR1Ca);
  IR1C = min(IR1Ca,2*IR2Da); IR2C = min(IR2Ca,2*IR1Da);
  
  % Blend into discharge voltage
  blend = (0:length(disV)-1)/(length(disV)-1);
  IRblend1 = IR1D + (IR2D-IR1D)*blend(:);
  Vdis = disV + IRblend1;
  
  % Blend into charge voltage
  blend = (0:length(chgV)-1)/(length(chgV)-1);
  IRblend3 = IR1C + (IR2C-IR1C)*blend(:);
  Vchg = chgV - IRblend3;
  
  deltaV50 = linearinterp(chgZ,Vchg,0.5) - linearinterp(disZ,Vdis,0.5);
  ind  = find(chgZ < 0.5);
  vChg = Vchg(ind) - chgZ(ind)*deltaV50;
  zChg = chgZ(ind);
  ind  = find(disZ > 0.5);
  vDis = flipud(Vdis(ind) + (1 - disZ(ind))*deltaV50);
  zDis = flipud(disZ(ind));
  
  Z  = stdZ;
  U  = linearinterp([zChg; zDis],[vChg; vDis],stdZ);  
end
  
function [Z, U] = OCP_voltageAverage(data, ~)
%OCP_VOLTAGEAVERAGE Average branch voltages at matching SOC values.
%   [Z, U] = OCP_voltageAverage(data, config) uses data.interp.stdZ and the
%   interpolated chgV/disV vectors. config is unused. Z is dimensionless;
%   U is the averaged voltage [V]. Both outputs are column vectors.
  
  Z = data.interp.stdZ;
  U = (data.interp.chgV+data.interp.disV)/2;
  
end
  
function [Z, U] = OCP_socAverage(data, ~)
%OCP_SOCAVERAGE Average branch SOC values at matching voltages.
%   [Z, U] = OCP_socAverage(data, config) averages data.interp.chgZ/disZ on
%   stdV, then interpolates to stdZ. config is unused. Z is dimensionless;
%   U is voltage [V]. Both outputs are column vectors.
  
  stdZ = data.interp.stdZ;
  stdV = data.interp.stdV;
  
  % Average soc under the same voltage
  z3 = (data.interp.chgZ+data.interp.disZ)/2;
  Z  = stdZ;
  U  = linearinterp(z3,stdV,stdZ);
  
end
  
function [Z, U] = OCP_diagonalAverage(data, config)
%OCP_DIAGONALAVERAGE Align branches using differential-capacity correlations.
%   [Z, U] = OCP_diagonalAverage(data, config) uses data.orig/data.interp and
%   config.daxcorrfiltv/daxcorrfiltz to estimate voltage and SOC shifts.
%   config.datype selects 'useDis', 'useChg', or 'useAvg' after alignment.
%   Z is the dimensionless SOC column grid; U is the estimated voltage [V].
  
  stdV = data.interp.stdV;
  disZ = data.interp.disZ;
  chgZ = data.interp.chgZ;
  
  stdZ = data.interp.stdZ;
  disV = data.interp.disV;
  chgV = data.interp.chgV;
  
  % First, find cross-correlations for U vs dZdV
  % That is, differential capacities versus U
  du = mean(diff(stdV));
  intervals = config.daxcorrfiltv(stdV);
  if ~iscell(intervals), intervals = {intervals}; end
  lagU = zeros(size(intervals));
  for k = 1:length(intervals)
    ind = intervals{k};
    dzdv1 = roughdiff(disZ(ind),stdV(ind));
    dzdv3 = roughdiff(chgZ(ind),stdV(ind));
    [c,lag] = xcorr(dzdv1,dzdv3); 
    lagU(k) = abs(lag(c==max(c))*du);
  end
  lagU = mean(lagU);
  
  % Second, find cross-correlations for Z vs dZdV (Only for the defined SOC region to focus)
  dz = mean(diff(stdZ));
  intervals = config.daxcorrfiltz(stdZ);
  if ~iscell(intervals), intervals = {intervals}; end
  lagZ = zeros(size(intervals));
  for k = 1:length(intervals)
    ind = intervals{k};
    dzdv1 = roughdiff(stdZ(ind),disV(ind));
    dzdv3 = roughdiff(stdZ(ind),chgV(ind));
    [c,lag] = xcorr(dzdv1,dzdv3); 
    lagZ(k) = abs(lag(c==max(c))*dz);
  end
  lagZ = mean(lagZ);
  
  % Shift both branches toward the equilibrium curve before interpolation.
  uuD = data.orig.disV+lagU/2;
  zzD = data.orig.disZ+lagZ/2;
  U4D = linearinterp(zzD,uuD,stdZ);
  uuC = data.orig.chgV-lagU/2;
  zzC = data.orig.chgZ-lagZ/2;
  U4C = linearinterp(zzC,uuC,stdZ);
  switch config.datype
    case 'useDis' % base off of discharge data
      U4 = U4D;
    case 'useChg' % base off of charge data
      U4 = U4C;
    case 'useAvg' % average the basis
      U4 = (U4C + U4D)/2;
  end
  
  Z = stdZ;
  U = U4;
end
