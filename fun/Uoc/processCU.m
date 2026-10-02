function cal = processCU(test_data, cfg, voltChannels, voltChannelInd)

chCap  = test_data.ch_capacity;
disCap = abs(test_data.dch_capacity);

scripts = struct();                     % loop to insert check-up data into struct
for k = 1:numel(cfg)
    ind = find(test_data.step == cfg(k).step ...
               & ...
               test_data.TC_Counter1 == cfg(k).counter);
    
    scripts(k).name  = cfg(k).name;
    scripts(k).ind   = ind;
    
    scripts(k).disAh = getAhSafe(disCap, ind); % reset capacities to index
    scripts(k).chgAh = getAhSafe(chCap,  ind);
end

% coulombic efficiency
totDisAh = sum(arrayfun(@(s) s.disAh(end), scripts));
totChgAh = sum(arrayfun(@(s) s.chgAh(end), scripts));
eta25    = min([totDisAh / totChgAh, 1]);
if eta25 == 0; eta25 = 1; end

for k = 1:numel(scripts)
  scripts(k).chgAh = scripts(k).chgAh.*eta25; % apply to charge data 
end

% capacity
Q25 = scripts(1).disAh(end) + scripts(2).disAh(end) - scripts(1).chgAh(end) - scripts(2).chgAh(end);

switch voltChannels{voltChannelInd}
  case 'MainVolt'
      minV = round(min(test_data.v)*100)/100;
      maxV = round(max(test_data.v)*100)/100;
      disV = test_data.v(scripts(1).ind);
      chgV = test_data.v(scripts(3).ind);
  case 'AuxVolt'
      minV = round(min(test_data.Aux_V1)*100)/100;
      maxV = round(max(test_data.Aux_V1)*100)/100;
      disV = test_data.Aux_V1(scripts(1).ind);
      chgV = test_data.Aux_V1(scripts(3).ind);
  case 'AuxMinusMain'
      minV = round(min(test_data.Aux_V1-test_data.v)*100)/100;
      maxV = round(max(test_data.Aux_V1-test_data.v)*100)/100;
      disV = test_data.Aux_V1(scripts(1).ind)-test_data.v(scripts(1).ind);
      chgV = test_data.Aux_V1(scripts(3).ind)-test_data.v(scripts(3).ind);
  case 'AuxPlusMain'
      minV = round(min(test_data.Aux_V1+test_data.v)*100)/100;
      maxV = round(max(test_data.Aux_V1+test_data.v)*100)/100;
      disV = test_data.Aux_V1(scripts(1).ind)+test_data.v(scripts(1).ind);
      chgV = test_data.Aux_V1(scripts(3).ind)+test_data.v(scripts(3).ind);
end

% DISCHARGE arrays
disZ = 1 - scripts(1).disAh/Q25;    % SOC of the CC phase
disI = test_data.i(scripts(1).ind);
disT = test_data.t(scripts(1).ind);

% CHARGE arrays
chgZ = scripts(3).chgAh/Q25;        % SOC of the CC phase
chgI = test_data.i(scripts(3).ind);
chgT = test_data.t(scripts(3).ind);

% obtain differential curves
[disZ, disV, disdZ] = smoothdiff(disZ, disV); % disV now vs SOC
[chgZ, chgV, chgdZ] = smoothdiff(chgZ, chgV); % chgV now vs SOC

% Ensure V is monotonicly increasing (not decreasing)
if disV(1) > disV(end)
disZ  = flip(disZ);
disV  = flip(disV);
disdZ = flip(disdZ);
end
if chgV(1) > chgV(end)
chgZ  = flip(chgZ);
chgV  = flip(chgV);
chgdZ = flip(chgdZ);
end

% Store processed data.
cal.TdegC         = 25;
cal.C_rate_actual = mean(abs([disI(:); chgI(:)]))/Q25;
cal.eta           = eta25;
cal.QAh           = Q25;
cal.disZ          = disZ;
cal.disV          = disV;
cal.disdZ         = disdZ;
cal.chgZ          = chgZ;
cal.chgV          = chgV;
cal.chgdZ         = chgdZ;
cal.disT          = disT - disT(1); % start from t=0
cal.disI          = disI;
cal.disV          = disV;
cal.chgT          = chgT - chgT(1); % start from t=0
cal.chgI          = chgI;
cal.chgV          = chgV;
cal.minV          = minV;
cal.maxV          = maxV;