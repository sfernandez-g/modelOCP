function test_data = getDataArbin_XLSX(file, step2incrementCounter, save_path)
% Optimized function to extract Arbin cycler data (.xlsx)
% Mondragon Unibertsitatea | 2024 (Updated 2026)

[filePath, fileBaseName, ~] = fileparts(file);
if nargin < 3 || isempty(save_path)
    matFile = fullfile(filePath, fileBaseName + ".mat");
else
    if ~isfolder(save_path), mkdir(save_path); end
    matFile = fullfile(save_path, fileBaseName + ".mat");
end

%% --- Check if MAT exists ---
if exist(matFile, 'file')
    loadedData = load(matFile, 'test_data');
    test_data = loadedData.test_data;
    fprintf('.mat archive found and loaded: %s\n', matFile);
    return;
end

%% --- Process Multiple Excel Files ---
baseName = regexprep(fileBaseName, '_Wb_\d+$', '');
excelFiles = dir(fullfile(filePath, baseName + '_Wb_*.xlsx'));
if isempty(excelFiles)
    % Fallback if the file doesn't follow the _Wb_N pattern
    fileNames = {fileBaseName + ".xlsx"};
else
    fileNames = sort({excelFiles.name});
end

test_data = struct();
test_data.raw = table();

%% --- Main Loop ---
for k = 1:length(fileNames)
    file_k = fullfile(filePath, fileNames{k});
    sheets = sheetnames(file_k);
    
    % Read data with preserved variable names
    raw_info = readtable(file_k, "Sheet", sheets(end-1), "VariableNamingRule", "preserve");
    raw_meas = readtable(file_k, "Sheet", sheets(end), "VariableNamingRule", "preserve");
    
    % Merge raw tables
    test_data.raw = [test_data.raw; raw_meas];
    
    % Extract Metadata only once (from the first file)
    if k == 1
        h = height(raw_info);
        test_data.Test_Name = raw_info.Properties.VariableNames{4};
        test_data.Start_Date = raw_info{h, 2};
        test_data.Channel    = raw_info{h, 1};
        test_data.Schedule   = raw_info{h, "Test Name"};
        test_data.Mass       = raw_info{h, "Var8"};
        test_data.Capacity_Nominal = raw_info{h, "Var10"};
    end
end

%% --- Vectorized Data Extraction ---
% Assigning to short variables for calculation
t    = test_data.raw.("Test Time (s)");
i    = test_data.raw.("Current (A)");
v    = test_data.raw.("Voltage (V)");
step = test_data.raw.("Step Index");

test_data.t = t;
test_data.i = i;
test_data.v = v;
test_data.step = step;

%% --- Dynamic Aux Fields (Temperature, Aux Voltages) ---
varNames = test_data.raw.Properties.VariableNames;
targetPrefixes = {"Aux_Temperature", "Aux_Voltage", "TC_Counter"};
shortPrefixes = {"T", "Aux_V", "TC_Counter"};

for p = 1:length(targetPrefixes)
    cols = varNames(startsWith(varNames, targetPrefixes{p}));
    for n = 1:length(cols)
        fieldName = sprintf("%s%d", shortPrefixes{p}, n);
        test_data.(fieldName) = test_data.raw.(cols{n});
    end
end

%% --- Vectorized Counter Correction ---
if nargin >= 2 && ~isempty(step2incrementCounter)
    % Increment counter only when entering the specified step
    step_changes = [0; diff(step)];
    test_data.TC_Counter1 = cumsum(step == step2incrementCounter & step_changes ~= 0);
end

%% --- Vectorized Capacity Calculations ---
dt = [0; diff(t)];          
Ah = (i .* dt) / 3600;      

test_data.capacity     = cumsum(Ah);
test_data.ch_capacity  = cumsum(Ah .* (i > 0));
test_data.dch_capacity = cumsum(abs(Ah .* (i < 0)));

% Step Capacity (Vectorized)
% Groups data by step change and applies cumsum to each group
step_ids = [0; cumsum(diff(step) ~= 0)];
test_data.step_capacity = cell2mat(accumarray(step_ids + 1, Ah, [], @(x) {cumsum(x)}));

%% --- Save Result ---
save(matFile, 'test_data');
fprintf('.mat archive generated and saved: %s\n', matFile);

end