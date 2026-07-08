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