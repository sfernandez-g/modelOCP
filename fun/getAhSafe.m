% GETAHSAFE Extract capacity data for a selected index range.
%
% -- Usage --
% Ah = getAhSafe(cap,ind)
%
% -- Input --
% cap = cumulative capacity vector [Ah]
% ind = indices defining the selected test segment
%
% -- Output --
% Ah = capacity vector for the selected segment, shifted so that
%      the first value is zero [Ah].
%
% If ind is empty, Ah is returned as 0.
%
% This function is used to safely extract charge or discharge capacity
% from a selected step of the test while resetting the capacity
% reference to the beginning of that step.
function Ah = getAhSafe(cap, ind)
    if isempty(ind)
        Ah = 0;
    else
        Ah = cap(ind);
        Ah = Ah - Ah(1);
    end