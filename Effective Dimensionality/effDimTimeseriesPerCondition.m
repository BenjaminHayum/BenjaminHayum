function [dataTable] = effDimTimeseriesPerCondition(measure,DataSet,doExclude, varargin)
% effDimAcrossSubjects Calculate effective dimensionality and plot for a
% specific condition
%
%   For optional arguments, see '+ecogutils/postHocOptionsParser.m' and the
%   help for +patientSummaries/stateAverage.m
%
%   There are some hardcoded state/color maps in this function that will
%   need to be modified depending on which data are to be plotted.
useEEGonly = false;

o = ecogutils.postHocOptionsParser(measure,DataSet,varargin{:});
if ~iscell(o.DataSet), o.DataSet = {o.DataSet}; end

ROI_table = ecogutils.loadROItable();

% doExclude = 7;
exclusionString = '';
if doExclude==7
    ROI_table = ROI_table(ROI_table.newROINum~=7,:); % Exclude PFC
    exclusionString = ' exPFC';
elseif doExclude==6
    ROI_table = ROI_table(ROI_table.newROINum~=6,:); % Exclude AudRel
    exclusionString = ' exAudRel';
elseif doExclude==9
    ROI_table = ROI_table(ROI_table.newROINum~=9,:); % Exclude Other
    exclusionString = ' exOther';
end

switch DataSet
    case 'emergence'
        colorMap = ...
            containers.Map({'U','emergence','S'},...
            {[139,0,0]./255,[0,128,0]./255,[0,0,255]./255});
    case 'sleep'

        colorMap = ...
            containers.Map({'WS','N1','N2','N3','R'},...
            {[139,0,0]./255,[0,128,0]./255,[0,0,255]./255,[128,128,128]./255,[255,165,0]./255});
    case {'proOR','dexOR'}
        colorMap = ...
            containers.Map({'WA','S','U'},...
            {[139,0,0]./255,[0,128,0]./255,[0,0,255]./255});
    case 'exercise'
        colorMap = ...
            containers.Map({'PreExercise','PostExercise'},...
            {[139,0,0]./255,[0,128,0]./255});
    case 'postIctal'
        colorMap = ...
            containers.Map({'RS','RSictal','RSictaldlrm'},...
            {[139,0,0]./255,[0,128,0]./255,[0,0,255]./255});
    otherwise
        error('Need to add this DataSet to state/color map list in effDimTimeseries.m')
end

%% Load data and map states
[dataTable] = loadGoodData(o);

[patientList,ia,ib] = unique(dataTable.patientID,'stable');

subplotsSize = ceil(sqrt(length(patientList)));

%% Eigenvalue sums/ratios by State 
% All Subjects' Eigenvalues in one 
effDims = struct();
allEmbeddings = struct();

for iPatient = 1:length(patientList)

    patientRows = ib==iPatient;

    if contains(o.measure,'wPLI_')
        freqs = fieldnames(dataTable.(o.measure){ia(iPatient)}{1}.wPLI_debias);
        if length(freqs)>1, error('Only 1 frequency supported.'), end
        extractData = @(r) cellfun(@(s) s.wPLI_debias.(freqs{1}),r,'UniformOutput',false);
    elseif contains(o.measure,'envCorrDBT')
        freqs = fieldnames(dataTable.(o.measure){ia(iPatient)}{1}.envCorr);
        if length(freqs)>1, error('Only 1 frequency supported.'), end
        extractData = @(r) cellfun(@(s) s.envCorr.(freqs{1}),r,'UniformOutput',false);
    else
        error('Unknown/unsupported connMeasure');
    end

    allSegs = cellfun(extractData,dataTable.(o.measure)(patientRows),'UniformOutput',false);
    allSegTimes = dataTable.segmentTimes(patientRows);

    useDim = 1 + ndims(allSegs{1}{1});
    fullindex = repmat({':'},useDim-1,1);
    allSegs = cat(1,allSegs{:});
    allSegs = cat(useDim,allSegs{:});
    allStates = cat(1,dataTable.states{patientRows});
    allSegTimes = cat(1,allSegTimes{:});

    if any(isnan(allSegTimes))
        allSegTimes(:) = (1:length(allSegTimes))-1000;
    else
        allSegTimes = (allSegTimes-allSegTimes(1))/60;
    end

    theseChannels = dataTable.ECoGchannels{ia(iPatient)};
    [sortIndex,newROICount,firstSort,secondSort] = ecogutils.chanSort(theseChannels,ROI_table);

    if useEEGonly
        warning('Using EEG rather than ECoG channels')
        sortIndex = find(strcmp({theseChannels.oldROI},'EEG'));
    end


    theseChannels = theseChannels(sortIndex);
    eigValRatios = nan(size(allStates));

    effDimen = nan(size(allStates));

    allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).adjMat = cell(length(allStates),1);
    allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).eVals = cell(length(allStates),1);
    allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).eVecs = cell(length(allStates),1);
    allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).states = allStates;
    allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).channels = theseChannels;
    allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).P = cell(length(allStates),1);

     for iSeg = 1:length(allStates)
        ephysAdj = squeeze(allSegs(sortIndex,sortIndex,iSeg));
        allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).adjMat{iSeg} = ephysAdj;


        %ephysAdj = exp(ephysAdj/0.5);

        ephysAdj(isnan(ephysAdj)) = 0;
        ephysAdj(ephysAdj<0) = 0;
        ephysAdj(logical(eye(length(ephysAdj)))) = 0;

        if ~any(isnan(ephysAdj),'all') && ~any(isinf(ephysAdj),'all') && ~all(ephysAdj==0,'all')

            try
            [eVals,eVecs,~,P] = DMEAnalysis.doDME(ephysAdj,floor(size(ephysAdj,1)/3),false);
            allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).eVals{iSeg} = eVals;
            allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).eVecs{iSeg} = eVecs;
            allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).P{iSeg} = P;
            %         [eVals,~,~,~] = DMEAnalysis.doDME(ephysAdj,floor(size(ephysAdj,1)-1),false);
            %              [eVals,~,~,~] = DMEAnalysis.doDME(ephysAdj,floor(size(ephysAdj,1)/6),false);
            eigValRatios(iSeg) = sum(abs(eVals(2:4)))/sum(abs(eVals(2:end)));



            eValSum = sum(abs(eVals(2:end)));
            effDimen(iSeg) = (eValSum.^2/(sum(eVals(2:end).^2)))/length(eVals(2:end));
            catch
                disp('bad segment')
            end


        else
            disp('Warning: adj matrix has nan/inf');
        end
    end

    [stateList,~,stateIDs] = unique(allStates);

    effDims.(['patient' dataTable.patientID{ia(iPatient)}]).effDimen = effDimen;
    effDims.(['patient' dataTable.patientID{ia(iPatient)}]).stateList = stateList;
    effDims.(['patient' dataTable.patientID{ia(iPatient)}]).stateIDs = stateIDs;

end

%% Plotting Using Previously Acquired Data
%Comments are only describing differences from effDimTimeseries.m

%Each figure has all patients for a specific condition. So, we iterate 
%through conditions first and for each condition we iterate through all patients 
colorKeys = keys(colorMap);
for iState = 1:length(colorKeys)
    h = figure('Position',[65 158 1485 792],'Name',measure);
    figure(h)

    thisColor = colorMap(colorKeys{iState});
    currState = colorKeys{1,iState};

    for iPatient = 1:length(patientList) 
        subplot(subplotsSize, subplotsSize, iPatient);
        title(dataTable.patientID{ia(iPatient)})
        hold on;

        patientRows = ib==iPatient;
        
        allStates = cat(1,dataTable.states{patientRows});

        
        %Creates an logical array where 1 if the state is equal to
        %currState and 0 if not. To be used to index into allSegTimes and
        %effDimens
        segIndices = zeros(length(allStates), 1);
        for iList = 1:length(allStates)
            if strcmp(allStates{iList}, currState)
                segIndices(iList) = 1;
            end
        end
        segIndices = logical(segIndices);

        allSegTimes = dataTable.segmentTimes(patientRows);
        allSegTimes = cat(1,allSegTimes{:});
        if any(isnan(allSegTimes))
            allSegTimes(:) = (1:length(allSegTimes))-1000;
        else
            allSegTimes = (allSegTimes-allSegTimes(1))/60;
        end
        
        % Makes sure that the indicies aren't all 0 which means this
        % patient has none of the condition we want. In this case we'd plot
        % nothing and move onto the next patient
        if any(segIndices)
            plot(allSegTimes(segIndices),effDims.(['patient' dataTable.patientID{ia(iPatient)}]).effDimen(segIndices),'o','color',thisColor,'MarkerFaceColor',thisColor)
        end
        xlabel('Time (minutes)')

        %Edit this to whatever approximate limits are for your data
        ylim([0.35 0.7])
    end

    legend(colorKeys{iState},'AutoUpdate','off');

end
%% Output
outPath = [getLocalPath('ECoGoutput') 'Embeddings' filesep];
if ~exist(outPath, 'dir')
    mkdir(outPath);
    warning(['Created new directory ' outPath]);
end

if isempty(o.OptionSet)
    optionString = '';
else
    optionString = ['-' o.OptionSet];
end
fileOut = [outPath strcat(o.DataSet{:}) ' - ' o.measure optionString exclusionString '.mat'];
save(fileOut,'allEmbeddings','-v7.3');

fileOut = [outPath 'effDims' strcat(o.DataSet{:}) ' - ' o.measure optionString exclusionString '.mat'];
save(fileOut,'effDims','-v7.3');


end