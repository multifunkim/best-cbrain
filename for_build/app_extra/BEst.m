function BEst(varargin)

try
    %% Parsing scheme
    p = inputParser;
    
    addParameter(p,'inputData',defaultVal('inputData'),checkValFunc('inputData'))
    addParameter(p,'outputName',defaultVal('outputName'),checkValFunc('outputName'))
    
    addParameter(p,'memMethod',defaultVal('memMethod'),checkValFunc('memMethod'))
    addParameter(p,'sensorsType',defaultVal('sensorsType'),checkValFunc('sensorsType'))
    addParameter(p,'reconstructionWindow',defaultVal('reconstructionWindow'),checkValFunc('reconstructionWindow'))
    addParameter(p,'baselineWindow',defaultVal('baselineWindow'),checkValFunc('baselineWindow'))
    addParameter(p,'baseline',defaultVal('baseline'),checkValFunc('baseline'))
    addParameter(p,'normalization',defaultVal('normalization'),checkValFunc('normalization'))
    addParameter(p,'clusteringMethod',defaultVal('clusteringMethod'),checkValFunc('clusteringMethod'))
    addParameter(p,'mspWindow',defaultVal('mspWindow'),checkValFunc('mspWindow'))
    addParameter(p,'mspThresholdMethod',defaultVal('mspThresholdMethod'),checkValFunc('mspThresholdMethod'))
    addParameter(p,'mspThreshold',defaultVal('mspThreshold'),checkValFunc('mspThreshold'))
    addParameter(p,'neighborhoodOrder',defaultVal('neighborhoodOrder'),checkValFunc('neighborhoodOrder'))
    addParameter(p,'spatialSmoothing',defaultVal('spatialSmoothing'),checkValFunc('spatialSmoothing'))
    addParameter(p,'activeMeanInit',defaultVal('activeMeanInit'),checkValFunc('activeMeanInit'))
    addParameter(p,'activeProbaInit',defaultVal('activeProbaInit'),checkValFunc('activeProbaInit'))
    addParameter(p,'lambdaInit',defaultVal('lambdaInit'),checkValFunc('lambdaInit'))
    addParameter(p,'activeProbaThreshold',defaultVal('activeProbaThreshold'),checkValFunc('activeProbaThreshold'))
    addParameter(p,'activeVarCoef',defaultVal('activeVarCoef'),checkValFunc('activeVarCoef'))
    addParameter(p,'inactiveVarCoef',defaultVal('inactiveVarCoef'),checkValFunc('inactiveVarCoef'))
    addParameter(p,'noiseCovMethod',defaultVal('noiseCovMethod'),checkValFunc('noiseCovMethod'))
    addParameter(p,'optimMethod',defaultVal('optimMethod'),checkValFunc('optimMethod'))
    addParameter(p,'useParallel',defaultVal('useParallel'),checkValFunc('useParallel'))
    
    p.KeepUnmatched = true;
    p.PartialMatching = false;
    
    parse(p,varargin{:})
    
    if ~isempty(p.UsingDefaults)
        disp('*** Trying defaults values for: ')
        disp(p.UsingDefaults)
    end
    if ~isempty(fieldnames(p.Unmatched))
        disp('*** Extra inputs:')
        disp(p.Unmatched)
    end
    
    
    %% Creates the head model(s) and options structures to run MEM
    S = load(p.Results.inputData);
    OPTIONS = be_main(p.Results.memMethod);
    
    
    % Sensors type
    if strcmp(p.Results.sensorsType,'EEG')
        OPTIONS.mandatory.DataTypes = {'EEG'};
        channels = S.Options.eegChannels;
    elseif strcmp(p.Results.sensorsType,'MEG')
        OPTIONS.mandatory.DataTypes = {'MEG'};
        channels = S.Options.megChannels;
    elseif strcmp(p.Results.sensorsType,'EEG+MEG')
        OPTIONS.mandatory.DataTypes = {'EEG','MEG'};
        channels = sort([S.Options.eegChannels,S.Options.megChannels]);
    else
        error('Unknown sensors type: %s', p.Results.sensorsType)
    end
    
    if numel(intersect(OPTIONS.mandatory.DataTypes,...
            S.Options.HeadModelSensorsType)) ~= numel(OPTIONS.mandatory.DataTypes)
        for k = 1:numel(OPTIONS.mandatory.DataTypes)
            if ~any(strcmp(OPTIONS.mandatory.DataTypes{k},S.Options.HeadModelSensorsType))
                error('Missing the sensors type %s in the head model',OPTIONS.mandatory.DataTypes{k})
            end
        end
    elseif numel(intersect(OPTIONS.mandatory.DataTypes,...
            unique(S.Options.ChannelsType))) ~= numel(OPTIONS.mandatory.DataTypes)
        for k = 1:numel(OPTIONS.mandatory.DataTypes)
            if ~any(strcmp(OPTIONS.mandatory.DataTypes{k},unique(S.Options.ChannelsType)))
                error('Missing the sensors type %s in the recording',OPTIONS.mandatory.DataTypes{k})
            end
        end
    end
    
    
    % Head model
    if strcmp(p.Results.sensorsType,'EEG')
        BE_HEADMODEL.Gain.matrix = S.Options.Gain(S.Options.eegChannels,:);
        BE_HEADMODEL.Gain.modality = 'EEG';
    elseif strcmp(p.Results.sensorsType,'MEG')
        BE_HEADMODEL.Gain.matrix = S.Options.Gain(S.Options.megChannels,:);
        BE_HEADMODEL.Gain.modality = 'MEG';
    elseif strcmp(p.Results.sensorsType,'EEG+MEG')
        BE_HEADMODEL.Gain(1).matrix = S.Options.Gain(S.Options.eegChannels,:);
        BE_HEADMODEL.Gain(1).modality = 'EEG';
        BE_HEADMODEL.Gain(2).matrix = S.Options.Gain(S.Options.megChannels,:);
        BE_HEADMODEL.Gain(2).modality = 'MEG';
    end
    BE_HEADMODEL.vertex_connectivity = S.Options.VertConn;
    
    
    % Data
    OPTIONS.mandatory.Data = S.Options.Recording(channels,:);
    OPTIONS.mandatory.ChannelTypes = S.Options.ChannelsType(channels);
    OPTIONS.mandatory.DataTime = S.Time;
    OPTIONS.optional.TimeSegment = getTimeWindow(p.Results.reconstructionWindow);
    
    
    % Baseline
    if isempty(p.Results.baseline)
        OPTIONS.optional.Baseline = OPTIONS.mandatory.Data;
        OPTIONS.optional.BaselineTime = OPTIONS.mandatory.DataTime;
    else
        OPTIONS.optional.Baseline = getfield(getfield(...
            load(p.Results.baseline,'Options'),'Options'),'Recording');
        OPTIONS.optional.Baseline = OPTIONS.optional.Baseline(channels,:);
        OPTIONS.optional.BaselineTime = getfield(load(p.Results.baseline,'Time'),'Time');
    end
    OPTIONS.optional.BaselineSegment = getTimeWindow(p.Results.baselineWindow);
    
    
    % Other options
    OPTIONS.optional.normalization = p.Results.normalization;
    OPTIONS.clustering.clusters_type = p.Results.clusteringMethod;
    OPTIONS.clustering.MSP_window = str2double(p.Results.mspWindow);
    if strcmp(p.Results.mspThresholdMethod,'arbitrary')
        OPTIONS.clustering.MSP_scores_threshold = str2double(p.Results.mspThreshold);
    else
        OPTIONS.clustering.MSP_scores_threshold = 'fdr';
    end
    OPTIONS.clustering.neighborhood_order = str2double(p.Results.neighborhoodOrder);
    OPTIONS.solver.spatial_smoothing = str2double(p.Results.spatialSmoothing);
    OPTIONS.model.active_mean_method = str2double(p.Results.activeMeanInit(1));
    OPTIONS.model.alpha_method = str2double(p.Results.activeProbaInit(1));
    OPTIONS.model.initial_lambda = str2double(p.Results.lambdaInit(1));
    OPTIONS.model.alpha_threshold = str2double(p.Results.activeProbaThreshold);
    OPTIONS.solver.active_var_mult = str2double(p.Results.activeVarCoef);
    OPTIONS.solver.inactive_var_mult = str2double(p.Results.inactiveVarCoef);
    OPTIONS.solver.NoiseCov_method = str2double(p.Results.noiseCovMethod(1));
    OPTIONS.solver.Optim_method = p.Results.optimMethod;
    OPTIONS.solver.parallel_matlab = strcmp(p.Results.useParallel,'true');
    
    
    % Clean output
    S = rmfield(S,'Options');
    
    
    %% Display
    disp('*** List of options: ')
    
    fprintf('Pipeline:\n')
    fprintf(['\t',OPTIONS.mandatory.pipeline,'\n\n'])
    
    fprintf('Data definition:\n')
    disp(struct('DataTypes',OPTIONS.mandatory.DataTypes,...
        'TimeSegment',OPTIONS.optional.TimeSegment,...
        'BaselineSegment',OPTIONS.optional.BaselineSegment,...
        'Baseline',OPTIONS.optional.Baseline,...
        'normalization',OPTIONS.optional.normalization,...
        'verbose',OPTIONS.optional.verbose))
    
    fprintf('Clustering:\n')
    disp(OPTIONS.clustering)
    
    fprintf('Model:\n')
    disp(OPTIONS.model)
    
    fprintf('Solver:\n')
    disp(OPTIONS.solver)
    
    
    %% Run the chosen inverse algorithm
    [R,OPTIONS] = be_main(BE_HEADMODEL,OPTIONS);
    S.ImageGridAmp = R.ImageGridAmp;
    S.MEMoptions = R.MEMoptions;
    if strcmp(p.Results.sensorsType,'EEG')
        S.Comment = [OPTIONS.Comment(1:4),': EEG',OPTIONS.Comment(5:end)];
    elseif strcmp(p.Results.sensorsType,'MEG')
        S.Comment = [OPTIONS.Comment(1:4),': MEG',OPTIONS.Comment(5:end)];
    elseif strcmp(p.Results.sensorsType,'EEG+MEG')
        S.Comment = [OPTIONS.Comment(1:4),': EEG-MEG',OPTIONS.Comment(5:end)];
    end
    save(p.Results.outputName,'-struct','S','-v7.3')
catch e
    disp(getReport(e))
    exit(1)
end
exit(0)
end


function val = validVal(argName)
switch argName
    case 'memMethod'
        val = {'cMEM'};
    case 'sensorsType'
        val = {'EEG';'MEG';'EEG+MEG'};
    case 'normalization'
        val = {'adaptive';'fixed'};
    case 'clusteringMethod'
        val = {'blockwise';'static'};
    case 'mspThresholdMethod'
        val = {'arbitrary';'fdr'};
    case 'activeMeanInit'
        val = {'1: regular minimum norm';...
            '2: null hypothesis';...
            '3: MSP-regularized minimum norm';...
            '4: L-curve optimized Minimum Norm Estimate'};
    case 'activeProbaInit'
        val = {'1: mean MSP scores';...
            '2: max MSP scores';...
            '3: median MSP scores';...
            '4: equal to 0.5';...
            '5: equal to 1'};
    case 'lambdaInit'
        val = {'0: null hypothesis (vector of zeros)';...
            '1: random'};
    case 'noiseCovMethod'
        val = {'0: Identity matrix';...
            '1: Scalar matrix';...
            '2: Diagonal matrix';...
            '3: Full';...
            '4: Wavelet-based'};
    case 'optimMethod'
        val = {'fminunc','minfunc'};
    case 'useParallel'
        val = {'true','false'};
    otherwise
        error('No list of valid arguments for ''%s'', type help ''help %s'' for more info.',argName,mfilename)
end
end


function val = defaultVal(argName)
switch argName
    case 'inputData'
        val = '';
    case 'outputName'
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
    otherwise
        error('No default argument for ''%s'', type help ''help %s'' for more info.',argName,mfilename)
end
end


function f = checkValFunc(argName)
switch argName
    case 'inputData'
        f = @(x) validateattributes(x,{'char'},{'nonempty'});
    case 'outputName'
        f = @(x) validateattributes(x,{'char'},{'nonempty'});
    case 'memMethod'
        f = @(x) checkString(x,validVal(argName));
    case 'sensorsType'
        f = @(x) checkString(x,validVal(argName));
    case 'reconstructionWindow'
        f = @(x) checkTimeWindow(x);
    case 'baselineWindow'
        f = @(x) checkTimeWindow(x);
    case 'baseline'
        f = @(x) validateattributes(x,{'char'});
    case 'normalization'
        f = @(x) checkString(x,validVal(argName));
    case 'clusteringMethod'
        f = @(x) checkString(x,validVal(argName));
    case 'mspWindow'
        f = @(x) validateattributes(str2double(x),{'double'},{'positive','finite','nonnan'});
    case 'mspThresholdMethod'
        f = @(x) checkString(x,validVal(argName));
    case 'mspThreshold'
        f = @(x) validateattributes(str2double(x),{'double'},{'>=',0,'<=',1});
    case 'neighborhoodOrder'
        f = @(x) checkPositiveInteger(x);
    case 'spatialSmoothing'
        f = @(x) validateattributes(str2double(x),{'double'},{'>=',0,'<=',1});
    case 'activeMeanInit'
        f = @(x) checkString(x,validVal(argName));
    case 'activeProbaInit'
        f = @(x) checkString(x,validVal(argName));
    case 'lambdaInit'
        f = @(x) checkString(x,validVal(argName));
    case 'activeProbaThreshold'
        f = @(x) validateattributes(str2double(x),{'double'},{'>=',0,'<=',1});
    case 'activeVarCoef'
        f = @(x) validateattributes(str2double(x),{'double'},{'>=',0,'<=',1});
    case 'inactiveVarCoef'
        f = @(x) validateattributes(str2double(x),{'double'},{'>=',0,'<=',1});
    case 'noiseCovMethod'
        f = @(x) checkString(x,validVal(argName));
    case 'optimMethod'
        f = @(x) checkString(x,validVal(argName));
    case 'useParallel'
        f = @(x) checkString(x,validVal(argName));
    otherwise
        error('No default argument for ''%s'', type help ''help %s'' for more info.',argName,mfilename)
end
end


function TF = checkString(str,validStrings)
if ~any(strcmp(str,validStrings))
    tmpStr = [];
    for i = 1:length(validStrings)
        tmpStr = [tmpStr,validStrings{i},', '];%#ok
    end
    tmpStr = tmpStr(1:(end-2));
    error('The input did not match any of these values:\n\n%s',tmpStr)
else
    TF = true;
end
end


function TF = checkTimeWindow(str)
exp = '(?<timeBegin>-?\d+(\.\d+)?)[^0-9-]+(?<timeEnd>-?\d+(\.\d+)?).*';
matchedStr = regexp(str,exp,'names');
if isempty(matchedStr)
    error('The input did not match the expected format: ''floatValue1 floatValue2'' (where the decimal place is given by the dot ''.'').')
end
tw = [matchedStr.timeBegin,' ',matchedStr.timeEnd];
validateattributes(str2double(strsplit(tw,' ')),{'double'},{'numel',2,'increasing','finite','nonnan'});
TF = true;
end


function tw = getTimeWindow(str)
exp = '(?<timeBegin>-?\d+(\.\d+)?)[^0-9-]+(?<timeEnd>-?\d+(\.\d+)?).*';
matchedStr = regexp(str,exp,'names');
tw = [str2double(matchedStr.timeBegin),str2double(matchedStr.timeEnd)];
end


function TF = checkPositiveInteger(str,varargin)
matchedStr = regexp(str,'^\d+$','match');
if isempty(matchedStr)||~strcmp(str,matchedStr)
    error('The input did not match the expected format: ''positiveIntegerValue''.')
end
validateattributes(str2double(str),{'double'},{'finite','nonnan'});
if (~isempty(varargin)&&~ismember(str2double(str),varargin{1}(:)))
    tmpStr = [];
    for i = 1:length(varargin{1})
        tmpStr = [tmpStr,int2str(varargin{1}(i)),', '];%#ok
    end
    tmpStr = tmpStr(1:(end-2));
    error('The input did not match any of these values:\n\n%s',tmpStr)
else
    TF = true;
end
end