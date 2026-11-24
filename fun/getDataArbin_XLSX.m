% GETDATAARBIN_XLSX Import and process Arbin battery test data from XLSX files.
%
% -- Usage --
% test_data = getDataArbin_XLSX(file)
% test_data = getDataArbin_XLSX(file,step2incrementCounter)
% test_data = getDataArbin_XLSX(file,step2incrementCounter,save_path)
% test_data = getDataArbin_XLSX(file,step2incrementCounter,save_path,forceReprocess)
%
% -- Input --
% file                  = path to one Arbin XLSX file. If multiple files
%                         with the same base name and different _Wb_N
%                         suffixes are present, they are automatically
%                         detected and concatenated.
% step2incrementCounter = step number used to reconstruct TC_Counter1
%                         when the TC_Counter provided by Arbin is not
%                         reliable (optional).
% save_path             = directory where the processed MAT file is saved
%                         (optional; defaults to .cache/ beside the XLSX files).
% forceReprocess        = true to ignore any cached MAT file (default false).
%
% -- Output --
% test_data = structure containing the processed test data, including:
%             .raw          = raw measurement table
%             .step         = step number vector; measurement vectors are N-by-1
%             .t             = test time [s]
%             .i             = current [A]
%             .v             = voltage [V]
%             .capacity      = cumulative signed capacity [Ah]
%             .step_capacity = signed capacity relative to each step [Ah]
%             .ch_capacity  = cumulative charge capacity [Ah]
%             .dch_capacity = cumulative discharge capacity [Ah]
%             .TC_Counter1  = cycle counter from Arbin or reconstructed
%                             from step transitions when requested
%
% The function also detects auxiliary temperature, voltage and TC_Counter
% columns and stores them as fields of test_data (e.g. T1, Aux_V1,
% TC_Counter1).
%
% When step2incrementCounter is provided, TC_Counter1 is reconstructed
% from changes in the specified step rather than using the counter
% provided by the cycler.

function test_data = getDataArbin_XLSX(file, step2incrementCounter, save_path, forceReprocess)

if nargin < 2, step2incrementCounter = []; end
if ~isempty(step2incrementCounter)
    validateattributes(step2incrementCounter, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'nonnegative'});
end
if nargin < 4 || isempty(forceReprocess), forceReprocess = false; end
if islogical(forceReprocess)
    validateattributes(forceReprocess, {'logical'}, {'scalar'});
else
    validateattributes(forceReprocess, {'numeric'}, {'scalar','real','finite','>=',0,'<=',1});
end
forceReprocess = logical(forceReprocess);

[filePath, fileBaseName, ~] = fileparts(file);
if nargin < 3 || isempty(save_path)
    save_path = fullfile(filePath, '.cache');
end
matFile = fullfile(save_path, fileBaseName + ".mat");

%% --- Process Multiple Excel Files ---
baseName = regexprep(fileBaseName, '_Wb_\d+$', '');
excelFiles = dir(fullfile(filePath, [baseName '_Wb_*.xlsx']));
if isempty(excelFiles)
    fileNames = {fileBaseName + ".xlsx"};
else
    wbNumbers = zeros(numel(excelFiles), 1);
    for k = 1:numel(excelFiles)
        token = regexp(excelFiles(k).name, '_Wb_(\d+)\.xlsx$', 'tokens', 'once');
        wbNumbers(k) = str2double(token{1});
    end
    [~, order] = sort(wbNumbers);
    fileNames = {excelFiles(order).name};
end

% Reuse the cache only when all source files and processing options match.
sourceFiles = cellfun(@(name) fullfile(filePath, name), fileNames, 'UniformOutput', false);
if any(~cellfun(@isfile, sourceFiles))
    error('getDataArbin_XLSX:MissingSource', 'Source workbook not found: %s', file);
end
sourceInfo = cellfun(@(sourceFile) dir(sourceFile), sourceFiles);
sourceMeta = struct('name', sourceFiles(:), ...
                    'bytes', num2cell([sourceInfo.bytes]'), ...
                    'datenum', num2cell([sourceInfo.datenum]'));
if exist(matFile, 'file') && ~forceReprocess
    loadedData = load(matFile, 'test_data', 'cache_meta');
    if isfield(loadedData, 'test_data') && isfield(loadedData, 'cache_meta') ...
            && isstruct(loadedData.cache_meta) ...
            && all(isfield(loadedData.cache_meta, {'schemaVersion', 'step2incrementCounter', 'sources'})) ...
            && loadedData.cache_meta.schemaVersion == 1 && ...
            isequaln(loadedData.cache_meta.step2incrementCounter, step2incrementCounter) && ...
            isequaln(loadedData.cache_meta.sources, sourceMeta)
        test_data = loadedData.test_data;
        fprintf('.mat archive found and up to date: %s\n', matFile);
        return;
    end
    fprintf('.mat archive is stale; reprocessing: %s\n', matFile);
end

%% --- Main Loop ---
test_data = struct();
rawTables = cell(numel(fileNames), 1);
for k = 1:length(fileNames)
    file_k = fullfile(filePath, fileNames{k});
    sheets = sheetnames(file_k);
    
    % Read data with preserved variable names
    raw_info = readtable(file_k, "Sheet", sheets(end-1), "VariableNamingRule", "preserve");
    raw_meas = readtable(file_k, "Sheet", sheets(end), "VariableNamingRule", "preserve");
    
    % Merge raw tables
    rawTables{k} = raw_meas;
    
    % Extract Metadata only once (from the first file)
    if k == 1
        h = height(raw_info);
        test_data.Test_Name  = raw_info.Properties.VariableNames{4};
        test_data.Start_Date = raw_info{h, 2};
        test_data.Channel    = raw_info{h, 1};
        test_data.Schedule   = raw_info{h, "Test Name"};
        test_data.Mass       = raw_info{h, "Var8"};
        test_data.Capacity_Nominal = raw_info{h, "Var10"};
    end
end
test_data.raw = vertcat(rawTables{:});

%% --- Vectorized Data Extraction ---
% Assigning to short variables for calculation
t    = test_data.raw.("Test Time (s)");
i    = test_data.raw.("Current (A)");
v    = test_data.raw.("Voltage (V)");
step = test_data.raw.("Step Index");
validateattributes(t, {'numeric'}, {'column', 'real', 'finite', 'nonempty'});
validateattributes(i, {'numeric'}, {'column', 'real', 'finite', 'numel', numel(t)});
if any(diff(t) < 0)
    error('getDataArbin_XLSX:TimeReset', ...
        'Test time decreases across samples or workbook boundaries; resolve the time reference before integration.');
end

test_data.t    = t;
test_data.i    = i;
test_data.v    = v;
test_data.step = step;

%% --- Dynamic Aux Fields (Temperature, Aux Voltages) ---
varNames = test_data.raw.Properties.VariableNames;
targetPrefixes = {"Aux_Temperature", "Aux_Voltage", "TC_Counter"};
shortPrefixes = {"T", "Aux_V", "TC_Counter"};

for p = 1:length(targetPrefixes)
    cols = varNames(startsWith(varNames, targetPrefixes{p}));
    for n = 1:length(cols)
        fieldName             = sprintf("%s%d", shortPrefixes{p}, n);
        test_data.(fieldName) = test_data.raw.(cols{n});
    end
end

%% --- Vectorized Counter Correction ---
if nargin >= 2 && ~isempty(step2incrementCounter)
    % Increment counter only when entering the specified step
    step_changes          = [0; diff(step)];
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
step_ids                = [0; cumsum(diff(step) ~= 0)];
test_data.step_capacity = cell2mat(accumarray(step_ids + 1, Ah, [], @(x) {cumsum(x)}));

%% --- Save Result ---
if ~isfolder(save_path), mkdir(save_path); end
cache_meta.schemaVersion = 1;
cache_meta.step2incrementCounter = step2incrementCounter;
cache_meta.sources = sourceMeta;
save(matFile, 'test_data', 'cache_meta');
fprintf('.mat archive generated and saved: %s\n', matFile);

end
