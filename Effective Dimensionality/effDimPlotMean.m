function [dataTable] = effDimPlotMean(measure,DataSet,varargin)
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


figure(h)
hold on;
ylim([0, 1]) % Inclusive of all means in delirium dataset
%ylim([0,0.08]) % Inclusive of all stdevs in delirium dataset
xlim([0, 4])
title('Mean ED All Subjects')

xlabel('Conditions')
%ylabel('stdev ED')
ylabel('Mean ED')

colorKeys = keys(colorMap);
for iState = 1:length(colorKeys)
    thisColor = colorMap(colorKeys{iState});
    plot(-20,-20,'o','color',thisColor);
end
legend(colorKeys,'AutoUpdate','off');


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
        
    patientMeanEDs = zeros(1,length(analysisStates));
    
    %Plot points
    for iState = 1:length(analysisStates)
        if all(ismember(analysisStates{iState},stateList) == 0)
            continue;
        end
        
        thisColor = colorMap(analysisStates{iState});
        segIndices = stateIDs==iState;
        %stdevED = nanstd(effDimen(segIndices)); 
        meanED = mean(effDimen(segIndices),'omitnan');
        patientMeanEDs(iState) = meanED;
        %plot(iState,stdevED,'.','color',thisColor)
        plot(iState, meanED,'o','color',thisColor)

            
        effDimSD.(['patient' dataTable.patientID{ia(iPatient)}]).(analysisStates{iState}) = std(effDimen(segIndices));
    end
    
    %Plot connecting lines
    indicesED = ones(1, length(patientMeanEDs));
    for iState = 1:length(analysisStates)
        indicesED(patientMeanEDs == 0) = 0;
    end

    %Finds first and second indicies you want to plot out of the
    %patientMeanEDs which may have 0 in it
    skipPlot = 0;
    for iED = 1:length(patientMeanEDs)
        if indicesED(iED) == 1 
            firstIndex = iED;
            secondIndex = firstIndex + 1;
            if secondIndex <= length(patientMeanEDs)
                for nED = (iED+1):length(patientMeanEDs)
                    if indicesED(nED) == 1
                        secondIndex = nED;
                        skipPlot = 0;
                        break;
                    else
                        skipPlot = 1;
                        continue;
                    end
                end
            else
                skipPlot = 1;
                break;
            end
        else
            continue;
        end
        
        %Now  we have the first and second index to plot
        if skipPlot == 0
            plot([firstIndex secondIndex],[patientMeanEDs(firstIndex) patientMeanEDs(secondIndex)], 'Color', [0.5, 0.5, 0.5])
        end
    end
end

end