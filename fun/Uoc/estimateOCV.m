%ESTIMATEOCP Estimate OCV-vs-SOC curve from dis/charge curves using the
% four methods.
%
% -- Usage --
% ocprel = ESTIMATEOCV(elec,config,arg)
%
% -- Input --
% elec   = struct of electrode data from process[OCP|OCV].m
% config = struct of configuration data for electode (see fields below)
% arg    = struct of arguments passed to top-level fitOCP function
%
% `config` is a struct with the following fields:
%
%   .vmin     = minimum operating voltage [V]
%   .vmax     = maximum operating voltage [V]
%   .datype*  = type of diagonal averaging ('useChg','useDis', or 'useAvg')
%   .daxcorrfiltz* = function used to limit SOC interval over which xcorr 
%                    of dis/charge dZ-vs-Z curves are evalulated (see below)
%   .daxcorrfiltv* = function used to limit voltage interval over which xcorr 
%                    of dis/charge dZ-vs-V curves are evalulated (see below)
%
%   The .daxcorrfiltz and .daxcorrfiltv functions accept vectors
%   of the electrode SOC or voltage, respectively, and return a vector of
%   logical indicies describing where the cross-correlation should be
%   computed. To use multiple distinct intervals and average the xcorr over
%   the intervals, return a cell array with an entry for each interval.
%
% -- Output --
% ocvrel = relative OCV-vs-composition estimates @ each temp (struct array)
%          (see fields below)
%
% The struct array `ocvrel` has the following fields:
%
%   .elecname  : electrode name ('NEG', 'POS', or 'FULL')
%   .TdegC     : temperature [degC]
%   .estimates : struct array with the following fields:
%     .Z       : relative composition vector
%     .U       : potential vector [V]
%     .dZ      : differential capacity vector [V]
%     .name    : name of estimator
%     .elecname: name of electrode ('NEG', 'POS', or 'FULL') (same as above)
%     .TdegC   : temperature (for convienence, same as above) [degC]
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

function ocvrel = estimateOCV(elec, config, OCPResolution, OCPEstimators)

  % Iterate temperatures.
  clear ocvrel;
  for kt = length(elec.cal):-1:1 % !!! backwards to struct array doesn't grow
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
    for ke = length(OCPEstimators):-1:1  % !!! backwards so struct array doesn't grow
      estName = OCPEstimators{ke};

      % Generate OCP estimate.
      estFnName = ['OCP_' estName];  % name of the MATLAB fn to make estimate
      [Z,U] = feval(estFnName, data, config);  % call MATLAB fn

      % !!! BUGFIX 2024.26.01: retain data only from vmin to vmax
      Z = linearinterp(U, Z, data.interp.stdV);
      U = data.interp.stdV;

      % Calculate differential capacity
      dZ = roughdiff(Z,U);
      % [Z, U, dZ] = smoothdiff(Z,U);

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

% OCP_RESISTANCEBLEND Method #1: "resistance-blending" 
function [Z, U] = OCP_resistanceBlend(data, ~)
  
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
  
  deltaV50 = interp1(chgZ,Vchg,0.5) - interp1(disZ,Vdis,0.5);
  ind  = find(chgZ < 0.5);
  vChg = Vchg(ind) - chgZ(ind)*deltaV50;
  zChg = chgZ(ind);
  ind  = find(disZ > 0.5);
  vDis = flipud(Vdis(ind) + (1 - disZ(ind))*deltaV50);
  zDis = flipud(disZ(ind));
  
  Z  = stdZ;
  U  = interp1([zChg; zDis],[vChg; vDis],stdZ,'linear','extrap');  
end
  
% OCP_VOLTAGEAVERAGE Method #2: voltage averaging
function [Z, U] = OCP_voltageAverage(data, ~)
  
  Z = data.interp.stdZ;
  U = (data.interp.chgV+data.interp.disV)/2;
  
end
  
% OCP_SOCAVERAGE Method #3: SOC averaging
function [Z, U] = OCP_socAverage(data, ~)
  
  stdZ = data.interp.stdZ;
  stdV = data.interp.stdV;
  
  % Average soc under the same voltage
  z3 = (data.interp.chgZ+data.interp.disZ)/2;
  Z  = stdZ;
  U  = interp1(z3,stdV,stdZ,'linear','extrap');
  
end
  
% OCP_DIAGONALAVERAGE Method #4: Diagonal interpolation
function [Z, U] = OCP_diagonalAverage(data, config)
  
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
  
  % Interpolate ocv (using discharge data only)
  uuD = data.orig.disV+lagU/2;
  zzD = data.orig.disZ+lagZ/2;
  U4D = interp1(zzD,uuD,stdZ,'linear','extrap');
  uuC = data.orig.chgV-lagU/2;
  zzC = data.orig.chgZ-lagZ/2;
  U4C = interp1(zzC,uuC,stdZ,'linear','extrap');
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