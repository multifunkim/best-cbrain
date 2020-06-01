function BEst(varargin)

exitCode = 0;
try
    %% Parsing scheme
    p = inputParser;
    
    addParameter(p, 'inputData', defaultVal('inputData'), ...
        checkValFunc('inputData'))
    addParameter(p, 'outputDirName', defaultVal('outputDirName'), ...
        checkValFunc('outputDirName'))
    
    addParameter(p, 'memMethod', defaultVal('memMethod'), ...
        checkValFunc('memMethod'))
    addParameter(p, 'sensorsType', defaultVal('sensorsType'), ...
        checkValFunc('sensorsType'))
    addParameter(p, 'reconstructionWindow', ...
        defaultVal('reconstructionWindow'), ...
        checkValFunc('reconstructionWindow'))
    addParameter(p, 'baselineWindow', defaultVal('baselineWindow'), ...
        checkValFunc('baselineWindow'))
    addParameter(p, 'baseline', defaultVal('baseline'), ...
        checkValFunc('baseline'))
    addParameter(p, 'normalization', defaultVal('normalization'), ...
        checkValFunc('normalization'))
    addParameter(p, 'clusteringMethod', defaultVal('clusteringMethod'), ...
        checkValFunc('clusteringMethod'))
    addParameter(p, 'mspWindow', defaultVal('mspWindow'), ...
        checkValFunc('mspWindow'))
    addParameter(p, 'mspThresholdMethod', defaultVal('mspThresholdMethod'), ...
        checkValFunc('mspThresholdMethod'))
    addParameter(p, 'mspThreshold', defaultVal('mspThreshold'), ...
        checkValFunc('mspThreshold'))
    addParameter(p, 'neighborhoodOrder', defaultVal('neighborhoodOrder'), ...
        checkValFunc('neighborhoodOrder'))
    addParameter(p, 'spatialSmoothing', defaultVal('spatialSmoothing'), ...
        checkValFunc('spatialSmoothing'))
    addParameter(p, 'activeMeanInit', defaultVal('activeMeanInit'), ...
        checkValFunc('activeMeanInit'))
    addParameter(p, 'activeProbaInit', defaultVal('activeProbaInit'), ...
        checkValFunc('activeProbaInit'))
    addParameter(p, 'lambdaInit', defaultVal('lambdaInit'), ...
        checkValFunc('lambdaInit'))
    addParameter(p, 'activeProbaThreshold', ...
        defaultVal('activeProbaThreshold'), ...
        checkValFunc('activeProbaThreshold'))
    addParameter(p, 'activeVarCoef', defaultVal('activeVarCoef'), ...
        checkValFunc('activeVarCoef'))
    addParameter(p, 'inactiveVarCoef', defaultVal('inactiveVarCoef'), ...
        checkValFunc('inactiveVarCoef'))
    addParameter(p, 'noiseCovMethod', defaultVal('noiseCovMethod'), ...
        checkValFunc('noiseCovMethod'))
    addParameter(p, 'optimMethod', defaultVal('optimMethod'), ...
        checkValFunc('optimMethod'))
    addParameter(p, 'useParallel', defaultVal('useParallel'), ...
        checkValFunc('useParallel'))
    addParameter(p, 'maxWorkers', defaultVal('maxWorkers'), ...
        checkValFunc('maxWorkers'))
    
    p.KeepUnmatched = true;
    p.PartialMatching = false;
    
    parse(p, varargin{:})
    
    if ~isempty(p.UsingDefaults)
        disp('*** Trying defaults values for: ')
        disp(p.UsingDefaults)
    end
    if ~isempty(fieldnames(p.Unmatched))
        disp('*** Extra inputs: ')
        disp(p.Unmatched)
    end
    
    
    %% Run MEM
    tmpDirRoot = fullfile(p.Results.outputDirName, '.tmp');
    inputData = readInputData(p.Results.inputData, tmpDirRoot);
    for k = 1 : length(inputData)
        try
            runMEM(inputData{k}, p.Results)
        catch
            fprintf('A problem occured with the input: \n%s\n', inputData{k})
            exitCode = 1;
        end
    end
    removeDir(tmpDirRoot)
    
    
catch e
    disp(getReport(e))
    exit(1)
end


exit(exitCode)
end


%% ------------ Local functions ------------------------------------------------

function val = validVal(argName)
switch argName
    case 'inputData'
        val = {'.mat'; '.tgz'; '.tar.gz'; '.tar'};
    case 'memMethod'
        val = {'cMEM'};
    case 'sensorsType'
        val = {'EEG'; 'MEG'; 'EEG+MEG'};
    case 'normalization'
        val = {'adaptive'; 'fixed'};
    case 'clusteringMethod'
        val = {'blockwise'; 'static'};
    case 'mspThresholdMethod'
        val = {'arbitrary'; 'fdr'};
    case 'activeMeanInit'
        val = {'1: regular minimum norm'; ...
            '2: null hypothesis'; ...
            '3: MSP-regularized minimum norm'; ...
            '4: L-curve optimized Minimum Norm Estimate'};
    case 'activeProbaInit'
        val = {'1: mean MSP scores'; ...
            '2: max MSP scores'; ...
            '3: median MSP scores'; ...
            '4: equal to 0.5'; ...
            '5: equal to 1'};
    case 'lambdaInit'
        val = {'0: null hypothesis (vector of zeros)'; ...
            '1: random'};
    case 'noiseCovMethod'
        val = {'0: Identity matrix'; ...
            '1: Scalar matrix'; ...
            '2: Diagonal matrix'; ...
            '3: Full'; ...
            '4: Wavelet-based'};
    case 'optimMethod'
        val = {'fminunc', 'minfunc'};
    case 'useParallel'
        val = {'true', 'false'};
    otherwise
        error(['No list of valid arguments for ''%s'', ', ...
            'type help ''help %s'' for more info.'], argName, mfilename)
end
end


function val = defaultVal(argName)
switch argName
    case 'inputData'
        val = '';
    case 'outputDirName'
        val = '';
    case 'memMethod'
        val = 'cMEM';
    case 'sensorsType'
        val = '';
    case 'reconstructionWindow'
        val = '';
    case 'baselineWindow'
        val = '';
    case 'baseline'
        val = '';
    case 'normalization'
        val = 'adaptive';
    case 'clusteringMethod'
        val = 'static';
    case 'mspWindow'
        val = '10';
    case 'mspThresholdMethod'
        val = 'arbitrary';
    case 'mspThreshold'
        val = '0';
    case 'neighborhoodOrder'
        val = '4';
    case 'spatialSmoothing'
        val = '0.6';
    case 'activeMeanInit'
        val = '2: null hypothesis';
    case 'activeProbaInit'
        val = '3: median MSP scores';
    case 'lambdaInit'
        val = '1: random';
    case 'activeProbaThreshold'
        val = '0';
    case 'activeVarCoef'
        val = '0.05';
    case 'inactiveVarCoef'
        val = '0';
    case 'noiseCovMethod'
        val = '2: Diagonal matrix';
    case 'optimMethod'
        val = 'fminunc';
    case 'useParallel'
        val = 'true';
    case 'maxWorkers'
        val = '12';
    otherwise
        error(['No default argument for ''%s'', ', ...
            'type help ''help %s'' for more info.'], argName, mfilename)
end
end


function f = checkValFunc(argName)
switch argName
    case 'inputData'
        f = @(x) checkInputData(x, validVal(argName));
    case 'outputDirName'
        f = @(x) validateattributes(x, {'char'}, {'nonempty'});
    case 'memMethod'
        f = @(x) checkString(x, validVal(argName));
    case 'sensorsType'
        f = @(x) checkString(x, validVal(argName));
    case 'reconstructionWindow'
        f = @(x) checkTimeWindow(x);
    case 'baselineWindow'
        f = @(x) checkTimeWindow(x);
    case 'baseline'
        f = @(x) ischar(x);
    case 'normalization'
        f = @(x) checkString(x, validVal(argName));
    case 'clusteringMethod'
        f = @(x) checkString(x, validVal(argName));
    case 'mspWindow'
        f = @(x) validateattributes(str2double(x), {'double'}, ...
            {'positive', 'finite', 'nonnan'});
    case 'mspThresholdMethod'
        f = @(x) checkString(x, validVal(argName));
    case 'mspThreshold'
        f = @(x) validateattributes(str2double(x), ...
            {'double'}, {'>=', 0, '<=', 1});
    case 'neighborhoodOrder'
        f = @(x) checkPositiveInteger(x);
    case 'spatialSmoothing'
        f = @(x) validateattributes(str2double(x), ...
            {'double'}, {'>=', 0, '<=', 1});
    case 'activeMeanInit'
        f = @(x) checkString(x, validVal(argName));
    case 'activeProbaInit'
        f = @(x) checkString(x, validVal(argName));
    case 'lambdaInit'
        f = @(x) checkString(x, validVal(argName));
    case 'activeProbaThreshold'
        f = @(x) validateattributes(str2double(x), ...
            {'double'}, {'>=', 0, '<=', 1});
    case 'activeVarCoef'
        f = @(x) validateattributes(str2double(x), ...
            {'double'}, {'>=', 0, '<=', 1});
    case 'inactiveVarCoef'
        f = @(x) validateattributes(str2double(x), ...
            {'double'}, {'>=', 0, '<=', 1});
    case 'noiseCovMethod'
        f = @(x) checkString(x, validVal(argName));
    case 'optimMethod'
        f = @(x) checkString(x, validVal(argName));
    case 'useParallel'
        f = @(x) checkString(x, validVal(argName));
    case 'maxWorkers'
        f = @(x) validateattributes(str2double(x), ...
            {'double'}, {'positive', 'finite', 'nonnan'});
    otherwise
        error(['No default argument for ''%s'', ', ...
            'type help ''help %s'' for more info.'], argName, mfilename)
end
end


function TF = checkInputData(str, validTypes)
TF = any(cellfun(@(s) strcmp(str(end - length(s) + 1 : end), s), validTypes));
if TF
    return
elseif ~exist(str, 'dir')
    error(['The input data is not a directory or a file with format: %s, ', ...
        'type help ''help %s'' for more info.'], strjoin(validTypes), mfilename)
else
    TF = true;
end
end


function TF = checkString(str, validStrings)
if ~any(strcmp(str, validStrings))
    tmpStr = [];
    for i = 1 : length(validStrings)
        tmpStr = [tmpStr, validStrings{i}, ', ']; %#ok
    end
    tmpStr = tmpStr(1 : (end - 2));
    error('The input did not match any of these values: \n\n%s', tmpStr)
else
    TF = true;
end
end


function TF = checkTimeWindow(str)
exp = '(?<timeBegin>-?\d+(\.\d+)?)[^0-9-]+(?<timeEnd>-?\d+(\.\d+)?).*';
matchedStr = regexp(str, exp, 'names');
if isempty(matchedStr)
    error(['The input did not match the expected format: ', ...
        '''floatValue1 floatValue2'' (where the decimal place is given by', ...
        'the dot ''.'').'])
end
tw = [matchedStr.timeBegin, ' ', matchedStr.timeEnd];
validateattributes(str2double(strsplit(tw, ' ')), {'double'}, ...
    {'numel', 2, 'increasing', 'finite', 'nonnan'});
TF = true;
end


function tw = getTimeWindow(str)
exp = '(?<timeBegin>-?\d+(\.\d+)?)[^0-9-]+(?<timeEnd>-?\d+(\.\d+)?).*';
matchedStr = regexp(str, exp, 'names');
tw = [str2double(matchedStr.timeBegin), str2double(matchedStr.timeEnd)];
end


function TF = checkPositiveInteger(str, varargin)
matchedStr = regexp(str, '^\d+$', 'match');
if isempty(matchedStr)||~strcmp(str, matchedStr)
    error(['The input did not match the expected format: ', ...
        '''positiveIntegerValue''.'])
end
validateattributes(str2double(str), {'double'}, {'finite', 'nonnan'});
if (~isempty(varargin)&&~ismember(str2double(str), varargin{1}(:)))
    tmpStr = [];
    for i = 1 : length(varargin{1})
        tmpStr = [tmpStr, int2str(varargin{1}(i)), ', ']; %#ok
    end
    tmpStr = tmpStr(1 : (end - 2));
    error('The input did not match any of these values: \n\n%s', tmpStr)
else
    TF = true;
end
end


function inputData = readInputData(p, tmpDirRoot)
if strcmp(p(end - 3 : end), '.mat')
    inputData = {p};
elseif any(cellfun(@(s) strcmp(p(end - length(s) + 1 : end), s), ...
        {'.tgz', '.tar.gz', '.tar'}))
    tempDir = tempname(tmpDirRoot);
    createDir(tempDir)
    inputData = untar(p, tempDir);
    
    if iscell(inputData)
        discard = false(size(inputData));
        for k = 1 : length(inputData)
            if ~strcmp(inputData{k}(end - 3 : end), '.mat')
                discard(k) = true;
            end
        end
        inputData(discard) = [];
    elseif ~strcmp(inputData(end - 3 : end), '.mat')
        error('The input archive does not contain ''.mat'' files.')
    else
        inputData = {inputData};
    end
else
    inputData = [];
    files = dir(p);
    for k = 1 : length(files)
        if ~files(k).isdir
            inputData = [inputData, readInputData(fullfile(files(k).folder, ...
                files(k).name), tmpDirRoot)]; %#ok
        end
    end
end
end


function createDir(dirPath)
if ~exist(dirPath, 'dir')
    [status, msg] = mkdir(dirPath);
    if ~status
        error('- Failed to create the directory: \n  %s\n  %s\n', dirPath, msg)
    end
end
end


function removeDir(dirPath)
if exist(dirPath, 'dir')
    [status, msg] = rmdir(dirPath, 's');
    if ~status
        error('- Failed to remove the directory: \n  %s\n  %s\n', dirPath, msg)
    end
end
end


function S = runMEM(inputData, opt)
%% Creates the head model(s) and options structures to run MEM
S = load(inputData);
OPTIONS = be_main(opt.memMethod);


% Sensors type
if strcmp(opt.sensorsType, 'EEG')
    OPTIONS.mandatory.DataTypes = {'EEG'};
    channels = S.Options.eegChannels;
elseif strcmp(opt.sensorsType, 'MEG')
    OPTIONS.mandatory.DataTypes = {'MEG'};
    channels = S.Options.megChannels;
elseif strcmp(opt.sensorsType, 'EEG+MEG')
    OPTIONS.mandatory.DataTypes = {'EEG', 'MEG'};
    channels = sort([S.Options.eegChannels, S.Options.megChannels]);
else
    error('Unknown sensors type: %s', opt.sensorsType)
end

if numel(intersect(OPTIONS.mandatory.DataTypes, ...
        S.Options.HeadModelSensorsType)) ~= numel(OPTIONS.mandatory.DataTypes)
    for k = 1 : numel(OPTIONS.mandatory.DataTypes)
        if ~any(strcmp(OPTIONS.mandatory.DataTypes{k}, ...
                S.Options.HeadModelSensorsType))
            error('Missing the sensors type %s in the head model', ...
                OPTIONS.mandatory.DataTypes{k})
        end
    end
elseif numel(intersect(OPTIONS.mandatory.DataTypes, ...
        unique(S.Options.ChannelsType))) ~= numel(OPTIONS.mandatory.DataTypes)
    for k = 1 : numel(OPTIONS.mandatory.DataTypes)
        if ~any(strcmp(OPTIONS.mandatory.DataTypes{k}, ...
                unique(S.Options.ChannelsType)))
            error('Missing the sensors type %s in the recording', ...
                OPTIONS.mandatory.DataTypes{k})
        end
    end
end


% Head model
if strcmp(opt.sensorsType, 'EEG')
    BE_HEADMODEL.Gain.matrix = S.Options.Gain(S.Options.eegChannels, :);
    BE_HEADMODEL.Gain.modality = 'EEG';
elseif strcmp(opt.sensorsType, 'MEG')
    BE_HEADMODEL.Gain.matrix = S.Options.Gain(S.Options.megChannels, :);
    BE_HEADMODEL.Gain.modality = 'MEG';
elseif strcmp(opt.sensorsType, 'EEG+MEG')
    BE_HEADMODEL.Gain(1).matrix = S.Options.Gain(S.Options.eegChannels, :);
    BE_HEADMODEL.Gain(1).modality = 'EEG';
    BE_HEADMODEL.Gain(2).matrix = S.Options.Gain(S.Options.megChannels, :);
    BE_HEADMODEL.Gain(2).modality = 'MEG';
end
BE_HEADMODEL.vertex_connectivity = S.Options.VertConn;


% Data
OPTIONS.mandatory.Data = S.Options.Recording(channels, :);
OPTIONS.mandatory.ChannelTypes = S.Options.ChannelsType(channels);
OPTIONS.mandatory.DataTime = S.Time;
OPTIONS.optional.TimeSegment = getTimeWindow(opt.reconstructionWindow);


% Baseline
if isempty(opt.baseline)
    OPTIONS.optional.Baseline = OPTIONS.mandatory.Data;
    OPTIONS.optional.BaselineTime = OPTIONS.mandatory.DataTime;
else
    OPTIONS.optional.Baseline = getfield(getfield(...
        load(opt.baseline, 'Options'), 'Options'), 'Recording');
    OPTIONS.optional.Baseline = OPTIONS.optional.Baseline(channels, :);
    OPTIONS.optional.BaselineTime = getfield(load(opt.baseline, 'Time'), ...
        'Time');
end
OPTIONS.optional.BaselineSegment = getTimeWindow(opt.baselineWindow);


% Other options
OPTIONS.optional.normalization = opt.normalization;
OPTIONS.clustering.clusters_type = opt.clusteringMethod;
OPTIONS.clustering.MSP_window = str2double(opt.mspWindow);
if strcmp(opt.mspThresholdMethod, 'arbitrary')
    OPTIONS.clustering.MSP_scores_threshold = str2double(opt.mspThreshold);
else
    OPTIONS.clustering.MSP_scores_threshold = 'fdr';
end
OPTIONS.clustering.neighborhood_order = str2double(opt.neighborhoodOrder);
OPTIONS.solver.spatial_smoothing = str2double(opt.spatialSmoothing);
OPTIONS.model.active_mean_method = str2double(opt.activeMeanInit(1));
OPTIONS.model.alpha_method = str2double(opt.activeProbaInit(1));
OPTIONS.model.initial_lambda = str2double(opt.lambdaInit(1));
OPTIONS.model.alpha_threshold = str2double(opt.activeProbaThreshold);
OPTIONS.solver.active_var_mult = str2double(opt.activeVarCoef);
OPTIONS.solver.inactive_var_mult = str2double(opt.inactiveVarCoef);
OPTIONS.solver.NoiseCov_method = str2double(opt.noiseCovMethod(1));
OPTIONS.solver.Optim_method = opt.optimMethod;
OPTIONS.solver.parallel_matlab = strcmp(opt.useParallel, 'true');


% Clean output
S = rmfield(S, 'Options');


%% Display
disp('*** List of options: ')

fprintf('Pipeline: \n')
fprintf(['\t', OPTIONS.mandatory.pipeline, '\n\n'])

fprintf('Data definition: \n')
disp(struct('DataTypes', OPTIONS.mandatory.DataTypes, ...
    'TimeSegment', OPTIONS.optional.TimeSegment, ...
    'BaselineSegment', OPTIONS.optional.BaselineSegment, ...
    'Baseline', OPTIONS.optional.Baseline, ...
    'normalization', OPTIONS.optional.normalization, ...
    'verbose', OPTIONS.optional.verbose))

fprintf('Clustering: \n')
disp(OPTIONS.clustering)

fprintf('Model: \n')
disp(OPTIONS.model)

fprintf('Solver: \n')
disp(OPTIONS.solver)


%% Run the chosen inverse algorithm
if strcmp(opt.useParallel, 'true') && isempty(gcp('nocreate'))
    parpool('local', str2double(opt.maxWorkers));
end
[R, OPTIONS] = be_main(BE_HEADMODEL, OPTIONS);
S.ImageGridAmp = R.ImageGridAmp;
S.MEMoptions = R.MEMoptions;
if strcmp(opt.sensorsType, 'EEG')
    S.Comment = [OPTIONS.Comment(1 : 4), ': EEG', OPTIONS.Comment(5 : end)];
elseif strcmp(opt.sensorsType, 'MEG')
    S.Comment = [OPTIONS.Comment(1 : 4), ': MEG', OPTIONS.Comment(5 : end)];
elseif strcmp(opt.sensorsType, 'EEG+MEG')
    S.Comment = [OPTIONS.Comment(1 : 4), ': EEG-MEG', OPTIONS.Comment(5 : end)];
end
[~, name] = fileparts(inputData);
save(fullfile(opt.outputDirName, sprintf('cbrain_%s_%s', name, ...
    datestr(datetime('now'), 'yymmddHHMMSSFFF'))), '-struct', 'S', '-v7.3')
end
