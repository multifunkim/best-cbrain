function BEst(varargin)
exitCode = 0;
try
    [inputData, outputDirName, memMethod, p] = getParam(varargin{:});
    createDir(outputDirName)
    tmpDirRoot = fullfile(outputDirName, '.tmp');
    inputData = readInputData(inputData, tmpDirRoot);
    for k = 1 : length(inputData)
        try
            runMEM(inputData{k}, outputDirName, memMethod, tmpDirRoot, ...
                p.Results)
        catch e
            fprintf('A problem occured when processing:\n%s\n%s\n', ...
                inputData{k}, getReport(e))
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

function [inputData, outputDirName, memMethod, p] = getParam(varargin)
% Primary arguments
p = inputParser;

addParameter(p, 'inputData', ...
    defaultVal('inputData', []), ...
    checkValFunc('inputData', []))
addParameter(p, 'outputDirName', ...
    defaultVal('outputDirName', []), ...
    checkValFunc('outputDirName', []))
addParameter(p, 'memMethod', ...
    defaultVal('memMethod'), ...
    checkValFunc('memMethod', []))

p.KeepUnmatched = true;
p.PartialMatching = false;

parse(p, varargin{:});

inputData = p.Results.inputData;
outputDirName = p.Results.outputDirName;
memMethod = p.Results.memMethod;


% Secondary arguments
p = inputParser;

addParameter(p, 'sensorsTypes', ...
    defaultVal('sensorsTypes', memMethod), ...
    checkValFunc('sensorsTypes', memMethod))
addParameter(p, 'reconstructionWindow', ...
    defaultVal('reconstructionWindow', memMethod), ...
    checkValFunc('reconstructionWindow', memMethod))
addParameter(p, 'baselineWindow', ...
    defaultVal('baselineWindow', memMethod), ...
    checkValFunc('baselineWindow', memMethod))
addParameter(p, 'baseline', ...
    defaultVal('baseline', memMethod), ...
    checkValFunc('baseline', memMethod))
addParameter(p, 'normalization', ...
    defaultVal('normalization', memMethod), ...
    checkValFunc('normalization', memMethod))
addParameter(p, 'clusteringMethod', ...
    defaultVal('clusteringMethod', memMethod), ...
    checkValFunc('clusteringMethod', memMethod))
addParameter(p, 'mspWindow', ...
    defaultVal('mspWindow', memMethod), ...
    checkValFunc('mspWindow', memMethod))
addParameter(p, 'mspThresholdMethod', ...
    defaultVal('mspThresholdMethod', memMethod), ...
    checkValFunc('mspThresholdMethod', memMethod))
addParameter(p, 'mspThreshold', ...
    defaultVal('mspThreshold', memMethod), ...
    checkValFunc('mspThreshold', memMethod))
addParameter(p, 'neighborhoodOrder', ...
    defaultVal('neighborhoodOrder', memMethod), ...
    checkValFunc('neighborhoodOrder', memMethod))
addParameter(p, 'spatialSmoothing', ...
    defaultVal('spatialSmoothing', memMethod), ...
    checkValFunc('spatialSmoothing', memMethod))
addParameter(p, 'activeMeanInit', ...
    defaultVal('activeMeanInit', memMethod), ...
    checkValFunc('activeMeanInit', memMethod))
addParameter(p, 'activeProbaInit', ...
    defaultVal('activeProbaInit', memMethod), ...
    checkValFunc('activeProbaInit', memMethod))
addParameter(p, 'lambdaInit', ...
    defaultVal('lambdaInit', memMethod), ...
    checkValFunc('lambdaInit', memMethod))
addParameter(p, 'activeProbaThreshold', ...
    defaultVal('activeProbaThreshold', memMethod), ...
    checkValFunc('activeProbaThreshold', memMethod))
addParameter(p, 'activeVarCoef', ...
    defaultVal('activeVarCoef', memMethod), ...
    checkValFunc('activeVarCoef', memMethod))
addParameter(p, 'inactiveVarCoef', ...
    defaultVal('inactiveVarCoef', memMethod), ...
    checkValFunc('inactiveVarCoef', memMethod))
addParameter(p, 'noiseCovMethod', ...
    defaultVal('noiseCovMethod', memMethod), ...
    checkValFunc('noiseCovMethod', memMethod))
addParameter(p, 'optimMethod', ...
    defaultVal('optimMethod', memMethod), ...
    checkValFunc('optimMethod', memMethod))
addParameter(p, 'frequencies', ...
    defaultVal('frequencies', memMethod), ...
    checkValFunc('frequencies', memMethod))
addParameter(p, 'waveletType', ...
    defaultVal('waveletType', memMethod), ...
    checkValFunc('waveletType', memMethod))
addParameter(p, 'vanishingMoments', ...
    defaultVal('vanishingMoments', memMethod), ...
    checkValFunc('vanishingMoments', memMethod))
addParameter(p, 'coefShrinkage', ...
    defaultVal('coefShrinkage', memMethod), ...
    checkValFunc('coefShrinkage', memMethod))
addParameter(p, 'useParallel', ...
    defaultVal('useParallel', memMethod), ...
    checkValFunc('useParallel', memMethod))
addParameter(p, 'maxWorkers', ...
    defaultVal('maxWorkers', memMethod), ...
    checkValFunc('maxWorkers', memMethod))

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
end


function val = validVal(argName, memMethod)
switch argName
    case 'inputData'
        val = {'.mat'; '.tgz'; '.tar.gz'; '.tar'; 'dir'};
        return
    case 'memMethod'
        val = {'cMEM'; 'wMEM'};
        return
    otherwise
        if isempty(memMethod)
            error(['No list of valid arguments for ''%s'', ', ...
                'MEM method not specified, ', ...
                'type help ''help %s'' for more info.'], argName, mfilename)
        end
end

switch argName
    case 'sensorsTypes'
        if strcmp(memMethod, 'cMEM')
            val = {'EEG'; 'MEG'; 'EEG+MEG'};
        elseif strcmp(memMethod, 'wMEM')
            val = {'EEG'; 'MEG'};
        end
    case 'baseline'
        val = {'.mat'; '.tgz'; '.tar.gz'; '.tar'};
    case 'normalization'
        val = {'adaptive'; 'fixed'};
        
    case 'clusteringMethod'
        if strcmp(memMethod, 'cMEM')
            val = {'blockwise'; 'static'};
        elseif strcmp(memMethod, 'wMEM')
            val = {'wavelet-adaptive'};
        end
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
        if strcmp(memMethod, 'cMEM')
            val = {'0: Identity matrix'; ...
                '1: Scalar matrix'; ...
                '2: Diagonal matrix'; ...
                '3: Full'; ...
                '4: Wavelet-based'};
        elseif strcmp(memMethod, 'wMEM')
            val = {'4: Diagonal matrix'; ...
                '5: Scalar matrix'};
        end
    case 'optimMethod'
        val = {'fminunc', 'minfunc'};
        
    case 'waveletType'
        if strcmp(memMethod, 'cMEM')
            val = '';
        elseif strcmp(memMethod, 'wMEM')
            val = {'rdw'};
        end
        
    case 'useParallel'
        val = {'true', 'false'};
        
    otherwise
        error(['No list of valid arguments for ''%s'', ', ...
            'type help ''help %s'' for more info.'], argName, mfilename)
end
end


function val = defaultVal(argName, memMethod)
switch argName
    case 'inputData'
        val = '';
        return
    case 'outputDirName'
        val = '';
        return
    case 'memMethod'
        val = '';
        return
    otherwise
        if isempty(memMethod)
            error(['No default argument for ''%s'', ', ...
                'MEM method not specified, ', ...
                'type help ''help %s'' for more info.'], argName, mfilename)
        end
end

switch argName
    case 'sensorsTypes'
        val = '';
    case 'reconstructionWindow'
        val = '';
    case 'baselineWindow'
        val = '';
    case 'baseline'
        val = '';
    case 'normalization'
        if strcmp(memMethod, 'cMEM')
            val = 'adaptive';
        elseif strcmp(memMethod, 'wMEM')
            val = 'fixed';
        end
        
    case 'clusteringMethod'
        if strcmp(memMethod, 'cMEM')
            val = 'static';
        elseif strcmp(memMethod, 'wMEM')
            val = 'wavelet-adaptive';
        end
    case 'mspWindow'
        if strcmp(memMethod, 'cMEM')
            val = '10';
        elseif strcmp(memMethod, 'wMEM')
            val = '';
        end
    case 'mspThresholdMethod'
        if strcmp(memMethod, 'cMEM')
            val = 'arbitrary';
        elseif strcmp(memMethod, 'wMEM')
            val = 'fdr';
        end
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
        if strcmp(memMethod, 'cMEM')
            val = '2: Diagonal matrix';
        elseif strcmp(memMethod, 'wMEM')
            val = '5: Scalar matrix';
        end
    case 'optimMethod'
        val = 'fminunc';
        
    case 'frequencies'
        if strcmp(memMethod, 'cMEM')
            val = '';
        elseif strcmp(memMethod, 'wMEM')
            val = 'all';
        end
        
    case 'waveletType'
        if strcmp(memMethod, 'cMEM')
            val = '';
        elseif strcmp(memMethod, 'wMEM')
            val = 'rdw';
        end
    case 'vanishingMoments'
        if strcmp(memMethod, 'cMEM')
            val = '';
        elseif strcmp(memMethod, 'wMEM')
            val = '4';
        end
    case 'coefShrinkage'
        if strcmp(memMethod, 'cMEM')
            val = '';
        elseif strcmp(memMethod, 'wMEM')
            val = '1';
        end
        
    case 'useParallel'
        val = 'true';
    case 'maxWorkers'
        val = '12';
        
    otherwise
        error(['No default argument for ''%s'', ', ...
            'type help ''help %s'' for more info.'], argName, mfilename)
end
end


function f = checkValFunc(argName, memMethod)
switch argName
    case 'inputData'
        f = @(x) checkInputData(x, validVal(argName, memMethod));
    case 'outputDirName'
        f = @(x) validateattributes(x, {'char'}, {'nonempty'});
    case 'memMethod'
        f = @(x) checkString(x, validVal(argName, memMethod));
        
    case 'sensorsTypes'
        f = @(x) checkString(x, validVal(argName, memMethod));
    case 'reconstructionWindow'
        f = @(x) checkTimeWindow(x);
    case 'baselineWindow'
        f = @(x) checkTimeWindow(x);
    case 'baseline'
        f = @(x) checkInputData(x, validVal(argName, memMethod));
    case 'normalization'
        f = @(x) checkString(x, validVal(argName, memMethod));
        
    case 'clusteringMethod'
        f = @(x) checkString(x, validVal(argName, memMethod));
    case 'mspWindow'
        f = @(x) validateattributes(str2double(x), {'double'}, ...
            {'positive', 'finite', 'nonnan'});
    case 'mspThresholdMethod'
        f = @(x) checkString(x, validVal(argName, memMethod));
    case 'mspThreshold'
        f = @(x) validateattributes(str2double(x), ...
            {'double'}, {'>=', 0, '<=', 1});
    case 'neighborhoodOrder'
        f = @(x) checkPositiveInteger(x, 0, []);
    case 'spatialSmoothing'
        f = @(x) validateattributes(str2double(x), ...
            {'double'}, {'>=', 0, '<=', 1});
        
    case 'activeMeanInit'
        f = @(x) checkString(x, validVal(argName, memMethod));
    case 'activeProbaInit'
        f = @(x) checkString(x, validVal(argName, memMethod));
    case 'lambdaInit'
        f = @(x) checkString(x, validVal(argName, memMethod));
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
        f = @(x) checkString(x, validVal(argName, memMethod));
    case 'optimMethod'
        f = @(x) checkString(x, validVal(argName, memMethod));
        
    case 'frequencies'
        f = @(x) checkFrequencies(x);
        
    case 'waveletType'
        f = @(x) checkString(x, validVal(argName, memMethod));
    case 'vanishingMoments'
        f = @(x) checkPositiveInteger(x, 0, 4);
    case 'coefShrinkage'
        f = @(x) validateattributes(str2double(x), ...
            {'double'}, {'>=', 0, '<=', 1});
        
    case 'useParallel'
        f = @(x) checkString(x, validVal(argName, memMethod));
    case 'maxWorkers'
        f = @(x) validateattributes(str2double(x), ...
            {'double'}, {'positive', 'finite', 'nonnan'});
        
    otherwise
        error(['No default argument for ''%s'', ', ...
            'type help ''help %s'' for more info.'], argName, mfilename)
end
end


function TF = checkInputData(str, validTypes)
TF = exist(str, 'file') && ...
    (any(strcmp('dir', validTypes)) || ...
    any(cellfun(@(s) strcmp(str(end - length(s) + 1 : end), s), validTypes)));

if ~TF
    error(['The input data is not an existing file with format: ''%s'', ', ...
        'type help ''help %s'' for more info.'], ...
        strjoin(validTypes, ' | '), mfilename)
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
expr = '(?<timeBegin>-?\d+(\.\d+)?)[^0-9-]+(?<timeEnd>-?\d+(\.\d+)?).*';
matchedStr = regexp(str, expr, 'names');
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
expr = '(?<timeBegin>-?\d+(\.\d+)?)[^0-9-]+(?<timeEnd>-?\d+(\.\d+)?).*';
matchedStr = regexp(str, expr, 'names');
tw = [str2double(matchedStr.timeBegin), str2double(matchedStr.timeEnd)];
end


function TF = checkPositiveInteger(str, S, L)
TF = all(isstrprop(str, 'digit'));
if ~TF
    error(['The input did not match the expected format: ', ...
        '''positiveIntegerValue''.'])
end

x = str2double(str);
if ~isempty(S) && (x < S)
    error('Expected input to be ''>= %d''.', S)
elseif ~isempty(L) && (x > L)
    error('Expected input to be ''<= %d''.', L)
end
end


function TF = checkFrequencies(str)
if isempty(str)
    error('No frequencies were specified.')
end

TF = true;
if strcmp('all', str)
    return
end

if ~isequal(1 : numel(str), regexp(str, '[0-9.;-]'))
    error(['When ''all'' is not used, ', ...
        'only digits and the following characters ''.;-'' are allowed'])
end

f = getFrequencies(str);
if isempty(f)
    error('Could not read frequencies: %s', str)
end
end


function f = getFrequencies(str)
f = [];
s = cellfun(@(s) strrep(s, '-', ';'), strsplit(str, ';'), ...
    'UniformOutput', false);
for k = 1 : numel(s)
    if isempty(s)
        continue
    end
    f_ = str2num(s{k}).'; %#ok
    if (numel(f_) == 1)
        f_ = [f_, f_]; %#ok
    end
    f = [f; f_]; %#ok
end
end


function scales = getWaveletScales(frequencies, t)
N = fix(log2(numel(t)));
if N < 7
    error('Not enough time samples')
end
N = N - min(N - 1, 3);

if strcmp(frequencies, 'all')
    scales = 1 : N;
else
    f = getFrequencies(frequencies);
    fb = (1 / (t(2) - t(1))) ./ (([2, 1].' * 2.^(1 : N)));
    
    if ((min(f(:)) < fb(1, end)) || (max(f(:)) > fb(2, 1)))
        error('Only frequencies in [%f, %f] can be analyzed.', ...
            fb(1, 2), fb(end, 1))
    end
    
    scales = [];
    for k = 1 : size(f, 1)
        s = find((fb(1, :) <= f(k, 1)) & (f(k, 1) <= fb(2, :)));
        if (f(k, 2) ~= f(k, 1))
            s = find((fb(1, :) <= f(k, 2)) & ...
                (f(k, 2) <= fb(2, :)), 1, 'last') : min(s);
        end
        scales = [scales, s]; %#ok
    end
    scales = unique(scales);
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


function runMEM(inputData, outputDirName, memMethod, tmpDirRoot, opt)
% Creates the head model(s) and options structures to run MEM
S = load(inputData);
OPTIONS = be_main(memMethod);


% Sensors types
if strcmp(opt.sensorsTypes, 'EEG')
    OPTIONS.mandatory.DataTypes = {'EEG'};
    channels = S.Options.eegChannels;
elseif strcmp(opt.sensorsTypes, 'MEG')
    OPTIONS.mandatory.DataTypes = {'MEG'};
    channels = S.Options.megChannels;
elseif strcmp(opt.sensorsTypes, 'EEG+MEG')
    OPTIONS.mandatory.DataTypes = {'EEG', 'MEG'};
    channels = sort([S.Options.eegChannels, S.Options.megChannels]);
else
    error('Unknown sensors types: %s', opt.sensorsTypes)
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
if strcmp(opt.sensorsTypes, 'EEG')
    BE_HEADMODEL.Gain.matrix = S.Options.Gain(S.Options.eegChannels, :);
    BE_HEADMODEL.Gain.modality = 'EEG';
elseif strcmp(opt.sensorsTypes, 'MEG')
    BE_HEADMODEL.Gain.matrix = S.Options.Gain(S.Options.megChannels, :);
    BE_HEADMODEL.Gain.modality = 'MEG';
elseif strcmp(opt.sensorsTypes, 'EEG+MEG')
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
    baseline = readInputData(opt.baseline, tmpDirRoot);
    if (numel(baseline) ~= 1)
        error('Expected a single baseline file but found ''%d'' in:\n%s', ...
            numel(baseline), opt.baseline)
    end
    OPTIONS.optional.Baseline = getfield(getfield(load(baseline{1}, ...
        'Options'), 'Options'), 'Recording');
    OPTIONS.optional.Baseline = OPTIONS.optional.Baseline(channels, :);
    OPTIONS.optional.BaselineTime = getfield(load(opt.baseline, 'Time'), ...
        'Time');
end
OPTIONS.optional.BaselineSegment = getTimeWindow(opt.baselineWindow);


% Other options
OPTIONS.optional.normalization = opt.normalization;
OPTIONS.clustering.clusters_type = opt.clusteringMethod;
if strcmp(memMethod, 'cMEM')
    OPTIONS.clustering.MSP_window = str2double(opt.mspWindow);
end
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
if strcmp(memMethod, 'wMEM')
    OPTIONS.wavelet.selected_scales = getWaveletScales(opt.frequencies, S.Time);
    OPTIONS.wavelet.type = opt.waveletType;
    OPTIONS.wavelet.vanish_moments = str2double(opt.vanishingMoments);
    OPTIONS.wavelet.shrinkage = str2double(opt.coefShrinkage);
end


% Clean output
S = rmfield(S, 'Options');


% Display
disp('*** List of options: ')

fprintf('Pipeline: \n')
fprintf(['\t', OPTIONS.mandatory.pipeline, '\n\n'])

fprintf('Data definition: \n')
disp(struct('DataTypes', {OPTIONS.mandatory.DataTypes}, ...
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

if strcmp(memMethod, 'wMEM')
    fprintf('Wavelet: \n')
    disp(OPTIONS.wavelet)
end


% Run the chosen inverse algorithm
if strcmp(opt.useParallel, 'true') && isempty(gcp('nocreate'))
    parpool('local', min(str2double(opt.maxWorkers), maxNumCompThreads));
end
[R, OPTIONS] = be_main(BE_HEADMODEL, OPTIONS);
S.ImageGridAmp = R.ImageGridAmp;
S.MEMoptions = R.MEMoptions;
if strcmp(opt.sensorsTypes, 'EEG')
    S.Comment = [OPTIONS.Comment(1 : 4), ': EEG', OPTIONS.Comment(5 : end)];
elseif strcmp(opt.sensorsTypes, 'MEG')
    S.Comment = [OPTIONS.Comment(1 : 4), ': MEG', OPTIONS.Comment(5 : end)];
elseif strcmp(opt.sensorsTypes, 'EEG+MEG')
    S.Comment = [OPTIONS.Comment(1 : 4), ': EEG-MEG', OPTIONS.Comment(5 : end)];
end
[~, name] = fileparts(inputData);
save(fullfile(outputDirName, sprintf('cbrain_%s_%s', name, ...
    datestr(datetime('now'), 'yymmddHHMMSSFFF'))), '-struct', 'S', '-v7.3')
end
