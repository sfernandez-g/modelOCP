function y = lithiationBalanceEquation(U, z, params, TdegC)
Z = compositionMSMR(U, params, TdegC);
y = z - Z;


