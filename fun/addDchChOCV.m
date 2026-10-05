% ADDDCHCHOCV Add experimental charge/discharge OCV data to an OCV structure.
%
% -- Usage --
% ocvrel = addDchChOCV(ocvrel,cellData,name)
%
% -- Input --
% ocvrel   = OCV structure containing the estimated OCV curves.
% cellData = structure containing the processed cell data from processCU.
% name     = curve to add:
%            "discharge" -> adds the discharge curve
%            "charge"    -> adds the charge curve
%
% -- Output --
% ocvrel = input OCV structure with the selected experimental curve
%          appended to ocvrel.estimates.
%
% The added estimate contains:
%          .Z       = composition/SOC vector [unitless]
%          .U       = voltage vector [V]
%          .dZ      = differential composition dZ/dU [1/V]
%          .name    = "discharge" or "charge"
%          .elecname = electrode/cell name
%          .TdegC   = test temperature [degC]

function ocvrel = addDchChOCV(ocvrel, cellData, name)

if name == "discharge"
    mini_name = "dis";
elseif name == "charge"
    mini_name = "chg";
end

ocvrel.estimates(end+1).Z = cellData.cal.(mini_name + "Z");
ocvrel.estimates(end).U = cellData.cal.(mini_name + "V");
ocvrel.estimates(end).dZ =cellData.cal.(mini_name + "dZ");
ocvrel.estimates(end).name = name;
ocvrel.estimates(end).elecname = cellData.name;
ocvrel.estimates(end).TdegC = cellData.cal.TdegC;