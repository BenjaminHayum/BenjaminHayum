function [dataTable] = effDimRSvsPostIctal(measure,DataSet,varargin)
%effDimTimeseries Calculate effective dimensionality and plot
%timeseries
%   Detailed explanation goes here

o = ecogutils.postHocOptionsParser(measure,DataSet,varargin{:});
if ~iscell(o.DataSet), o.DataSet = {o.DataSet}; end

ROI_table = ecogutils.loadROItable();

% ROI_table = ROI_table(ROI_table.newROINum~=7,:); % Exclude PFC
% ROI_table = ROI_table(ROI_table.newROINum~=6,:); % Exclude AudRel
% ROI_table = ROI_table(ROI_table.newROINum~=9,:); % Exclude Other


% colorMap = ...
%     containers.Map({'WS','N1','N2','N3','R'},...
%     {[139,0,0]./255,[0,128,0]./255,[0,0,255]./255,[128,128,128]./255,[255,165,0]./255});

% colorMap = ...
%     containers.Map({'WA','S','U'},...
%     {[139,0,0]./255,[0,128,0]./255,[0,0,255]./255});

% colorMap = ...
%    containers.Map({'PreExercise','PostExercise'},...
%    {[139,0,0]./255,[0,128,0]./255});

colorMap = ...
    containers.Map({'Ctrl','Ictal','IctalDlrm'},...
    {[139,0,0]./255,[0,128,0]./255,[0,0,255]./255});

% colorMap = ...
%     containers.Map({'S','U'},...
%     {[0,128,0]./255,[0,0,255]./255});

%% Load data and map states
[dataTable] = loadGoodData(o);

[patientList,ia,ib] = unique(dataTable.patientID,'stable');

%ia and ib are definitely the indexers for which comparisons are made!

h = figure('Position',[65 158 1485 792]);

%% Eigenvalue sums/ratios by state
effDims = struct();
effDimSD = struct();
allEmbeddings = struct();

colorKeys = keys(colorMap);

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
    
    useDim = 1 + ndims(allSegs{1}{1});
    fullindex = repmat({':'},useDim-1,1);
    allSegs = cat(1,allSegs{:});
    allSegs = cat(useDim,allSegs{:});
    allStates = cat(1,dataTable.states{patientRows});
    
    
    theseChannels = dataTable.ECoGchannels{ia(iPatient)};
    [sortIndex,newROICount,firstSort,secondSort] = ecogutils.chanSort(theseChannels,ROI_table);

    theseChannels = theseChannels(sortIndex);
    eigValRatios = nan(size(allStates));
    
    effDimen = nan(size(allStates));
    
    allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).adjMat = cell(length(allStates),1);
    allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).eVals = cell(length(allStates),1);
    allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).eVecs = cell(length(allStates),1);
    allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).states = allStates;
    allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).channels = theseChannels;
    
    for iSeg = 1:length(allStates)
        ephysAdj = squeeze(allSegs(sortIndex,sortIndex,iSeg));
        allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).adjMat{iSeg} = ephysAdj;
        
        
        %ephysAdj = exp(ephysAdj/0.5);
        
        ephysAdj(isnan(ephysAdj)) = 0;
        ephysAdj(ephysAdj<0) = 0;
        ephysAdj(logical(eye(length(ephysAdj)))) = 0;
        
        if iPatient == 9 && iSeg == 20
            break;
        end

        if ~any(isnan(ephysAdj),'all') && ~any(isinf(ephysAdj),'all') && ~all(ephysAdj==0,'all')
            [eVals,eVecs,~,~] = DMEAnalysis.doDME(ephysAdj,floor(size(ephysAdj,1)/3),false);
            allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).eVals{iSeg} = eVals;
            allEmbeddings.(['patient' dataTable.patientID{ia(iPatient)}]).eVecs{iSeg} = eVecs;
            %         [eVals,~,~,~] = DMEAnalysis.doDME(ephysAdj,floor(size(ephysAdj,1)-1),false);
            %              [eVals,~,~,~] = DMEAnalysis.doDME(ephysAdj,floor(size(ephysAdj,1)/6),false);
            eigValRatios(iSeg) = sum(abs(eVals(2:4)))/sum(abs(eVals(2:end)));
            
            eValSum = sum(abs(eVals(2:end)));
            %effDimen(iSeg) = prod((abs(eVals(2:end))./eValSum).^(-(abs(eVals(2:end))./eValSum)));
            
            
            
            effDimen(iSeg) = (eValSum.^2/(sum(eVals(2:end).^2)))/length(eVals(2:end));
            
            
            
        else
            disp('Warning: adj matrix has nan/inf');
        end
    end
    
    [stateList,~,stateIDs] = unique(allStates);

    effDims.(['patient' dataTable.patientID{ia(iPatient)}]).effDimen = effDimen;
    effDims.(['patient' dataTable.patientID{ia(iPatient)}]).stateList = stateList;
    effDims.(['patient' dataTable.patientID{ia(iPatient)}]).stateIDs = stateIDs; 

    analysisStates = keys(colorMap);

    %Convert stateIDs indicies to be for all conditions and not just those
    %on stateList
    for iState = 1:length(analysisStates)
        for iStateList = 1:length(stateList)
            if strcmp(analysisStates{iState}, stateList{iStateList})
                stateIDs(stateIDs == iStateList) = iState;
            end
        end
    end
    
    effDims.(['patient' dataTable.patientID{ia(iPatient)}]).effDimen = effDimen;
    effDims.(['patient' dataTable.patientID{ia(iPatient)}]).stateList = stateList;
    effDims.(['patient' dataTable.patientID{ia(iPatient)}]).stateIDs = stateIDs; 
        
    for iState = 1:length(analysisStates)
        segIndices = stateIDs==iState;
        effDimSD.(['patient' dataTable.patientID{ia(iPatient)}]).(analysisStates{iState}) = std(effDimen(segIndices));
    end
    
end

%Plotting
figure(h)
hold on;

ylim([0, 1])
xlim([0, 1])
title('Mean ED RS vs Ictal')
xlabel('Mean ED RS')
ylabel('Mean ED Ictal')

for iPatient = 1:length(patientList)
    effDimen = effDims.(['patient' dataTable.patientID{ia(iPatient)}]).effDimen;
    stateIDs = effDims.(['patient' dataTable.patientID{ia(iPatient)}]).stateIDs;

    meanCtrl = 0;
    meanIctal = 0;
    meanIctalDlrm = 0;

    analysisStates = keys(colorMap);
    for iState = 1:length(analysisStates)
        segIndices = stateIDs == iState;
        if iState == 1
            meanCtrl = mean(effDimen(segIndices),'omitnan');
        elseif iState == 2
            meanIctal = mean(effDimen(segIndices),'omitnan');
        else
            meanIctalDlrm = mean(effDimen(segIndices),'omitnan');
        end
    end

    if meanCtrl == 0
        continue;
    end
    if meanIctal ~= 0
        plot(meanCtrl, meanIctal, 'go')
    end
    if meanIctalDlrm ~= 0
        plot(meanCtrl, meanIctalDlrm, 'rx')
    end
end
plot([0,1] , [0,1] , 'b--', 'LineWidth', 1)
legend('Non-Delirious', 'Delirious')
pbaspect([1 1 1])
end