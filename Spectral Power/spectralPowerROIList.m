function [dataTable] = spectralPowerROIList(measure,DataSet,varargin)
% Finds the average spectral power for each ROI and sorts them across
% predefined conditions

%varargin is the classic postHocOptionsParser parameters until the last two
%which are 'Band' and the frequency band you are looking into

%% Load data table and ROI table

% Checking if 'Band' is an argument passed through and adjusting "o"
% correspondingly
bandBoolean = 0;
for i=1:length(varargin)
    if isequal('Band', varargin{i})
        bandBoolean = 1;
        break
    end
end
if bandBoolean == 1
    o = ecogutils.postHocOptionsParser(measure,DataSet,varargin{1:end-2});
else 
    o = ecogutils.postHocOptionsParser(measure,DataSet,varargin{:});
end
if ~iscell(o.DataSet), o.DataSet = {o.DataSet}; end


ROI_table = ecogutils.loadROItable();

% ROI_table = ROI_table(ROI_table.newROINum~=7,:); % Exclude PFC
% ROI_table = ROI_table(ROI_table.newROINum~=6,:); % Exclude AudRel
% ROI_table = ROI_table(ROI_table.newROINum~=9,:); % Exclude Other


[dataTable] = loadGoodData(o);

[patientList,ia,ib] = unique(dataTable.patientID,'stable');


%% Mapping Conditions and getting Band we want

%Find index of 'StateMap'
for i=1:length(varargin)
    if isequal('StateMap', varargin{i})
        stateMapIndex = i;
        break
    end
end

%Mapping the original conditions onto the new ones
dataTable.mapCondition = cell(length(dataTable.conditions), 1);
for i = 1:length(dataTable.conditions)
    dataTable.mapCondition{i} = varargin{stateMapIndex+1}(dataTable.conditions{i});
end

allConditions = unique(dataTable.mapCondition);

% For post-ictal:
%   RS vs ictal -- no matter whether ictal is delirious or non-delirious

% For post-operative:
%   post-op vs resolved OR post-op vs ctrl early OR ...

% Others...?


%Should always be the case...
if bandBoolean == 1
    band = varargin{end};
end

% Matching bands to frequencies
bandFreq = zeros(1,2);
if strcmp(band, 'Delta') == 1
    bandFreq = [0 4];
elseif strcmp(band, 'Theta') == 1
    bandFreq = [4 8];
elseif strcmp(band, 'Alpha') == 1
    bandFreq = [8 14];
elseif strcmp(band, 'Beta') == 1
    bandFreq = [14 30];
elseif strcmp(band, 'Gamma') == 1
    bandFreq = [30 50];
elseif strcmp(band, 'High Gamma') == 1
    bandFreq = [70 110];
elseif strcmp(band, 'All Gamma') == 1
    bandFreq = [30 120];
end

%% Adding each subject to ROI List

% Spectral powers to look through for ALL SUBJECTS
spectralCells = dataTable.specAnalysis;

% Make each patient/condition combo (i.e. block) have a dictionary of ROI
% key to average spectral power value
dictionaryCells = cell(length(dataTable.mapCondition), 2);

dictionaryCounter = 0;
for iPatient = 1:length(patientList)
    patientRows = ib==iPatient;
    patientSpecCells = spectralCells(patientRows);
    
    patientConditions = dataTable.mapCondition(patientRows);
    for iCondition = 1:length(patientConditions)
        conditionMap = containers.Map;
        dictionaryCounter = dictionaryCounter + 1;
        currPatient = patientList{iPatient};
        currCondition = patientConditions{iCondition};
        
        % To be used as the key in the dictionary/map
        currCombo = currPatient + " " + currCondition;
        
        %Get row indicies from ROIs
        currECoGchannels = dataTable.ECoGchannels{dictionaryCounter};
        ROIs = cell(1, length(currECoGchannels));
        for i = 1:length(currECoGchannels)
            ROIs{i} = currECoGchannels(i).oldROI;
        end

        %Now, for each ROI, get their indicies and index 
        allROIs = unique(ROIs);
        for iROI = 1:length(allROIs)
            %Now, index into spectralCells and append to allValues
            currROI = allROIs{iROI};
            rowIndices = [];
            for nROI = 1:length(ROIs)
                if strcmp(ROIs(nROI), currROI) == 1
                    rowIndices = [rowIndices nROI];
                end
            end
            
            allValues = [];
            currSpecCell = patientSpecCells{iCondition};
            specNames = fieldnames(currSpecCell);
            for iSegStruct = 1:size(specNames)
                currSegStruct = currSpecCell.(specNames{iSegStruct});
            
                % Based off of column indices in "freq" that correspond to the
                % frequencies recorded in "powspctrm"
                columnIndices = [];
                for i = 1:length(currSegStruct.freq)
                    if (currSegStruct.freq(i) >= bandFreq(1)) && (currSegStruct.freq(i) <= bandFreq(2))
                        columnIndices = [columnIndices i];
                    end
                end

                specValues = currSegStruct.powspctrm(rowIndices, columnIndices);
                %Turns it to be 1D
                specValues = specValues(:);
                allValues = [allValues specValues];
            end

            allValues = allValues(:);
            conditionMap(currROI) = mean(allValues);
            
        %end of allROIs loop
        end
        
        dictionaryCells{dictionaryCounter, 1} = currCombo;
        dictionaryCells{dictionaryCounter, 2} = conditionMap;

    %end of conditions loop
    end

%end of patients loop
end

%% Putting it all into a nicely formatted table

%Table for each patient + condition where in it is all the ROIs (keys) and
%average spectral powers (values)
spectralROIList = table();

%Below is for line 188
maxMapLength = 0;
for iDict = 1:length(dictionaryCells)
    mapLength = size(dictionaryCells{iDict, 2}, 1);
    if mapLength > maxMapLength
        maxMapLength = mapLength;
    end
end

%Sorting the ROIs by spectral power
for iPatient = 1:length(dictionaryCells)
    k = keys(dictionaryCells{iPatient, 2});
    v = cell2mat(values(dictionaryCells{iPatient, 2}));
    avg2roi = containers.Map(v, k);
    
    vSorted = sort(v, 'descend');
    v = num2cell(vSorted);
    
    spectralROIList.(dictionaryCells{iPatient, 1}) = cell(maxMapLength, 2);
    for iVal = 1:length(v)
        spectralROIList.(dictionaryCells{iPatient, 1}){iVal, 1} = avg2roi(v{iVal});
        spectralROIList.(dictionaryCells{iPatient, 1}){iVal, 2} = v{iVal};
    end
    
end

%% All Spectral Power Table Output
%
%fileOut = ['C:\Users\Benjamin Hayum\Desktop\spectralPowerROIList\' band 'PostIctalSpectralPowerROI.mat'];
%save(fileOut, 'spectralROIList');
%

%% Graphing

%If they are off of the diagonal then the ROI changes a lot and we should look
%into it

%Do one set of subplots of x = control spectral power and y = ictal spectral power
%with color of y changing according to delirious vs not delirious

h = figure('Position',[65 158 1485 792]);
subplotsSize = ceil(sqrt(length(patientList)));

%Horrible logic for this first set of graphs...
newPatientBoolean = 0;
patientCounter = 0;
currPatient = "ain't nobody";
for iCombo = 1:(length(dictionaryCells) + 1)
    lastRunBoolean = 0;
    if iCombo == (length(dictionaryCells) + 1)
        lastRunBoolean = 1;
        newPatientBoolean = 1;
    else
        newPatient = extractBetween(dictionaryCells{iCombo, 1}, 1, 4);
    end
    if strcmp(newPatient, currPatient) == 0 || lastRunBoolean == 1
        if patientCounter ~= 0
            newPatientBoolean = 1;
        end
        % Plot the old patient before we transition to the new 
        % (as long as the old patient exists which it won't if this is the 
        % first rum through --> why the below if statement is there)
        if newPatientBoolean == 1
            figure(h)
            subplot(subplotsSize,subplotsSize, patientCounter);
            hold on;
            % Error is "invalid parameter/value pair arguments"
            if ictalBoolean == 1
                text(ctrlSpectral, ictalSpectral, currROIsArray, 'Color', 'red', 'FontSize', 7, 'HorizontalAlignment', 'center', 'DisplayName', 'Non-Delirious');
                ictalMin = min(ictalSpectral);
                ictalMax = max(ictalSpectral);
            end
            if dlrmBoolean == 1
                text(ctrlSpectral, dlrmSpectral, currROIsArray, 'Color', 'green', 'FontSize', 7,'HorizontalAlignment', 'center', 'DisplayName', 'Delirious');
                dlrmMin = min(dlrmSpectral);
                dlrmMax = max(dlrmSpectral);
            end

            %Change the limits to the minimum and maximums
            currMin = min(ctrlSpectral);
            currMax = max(ctrlSpectral);
            if ictalBoolean == 1 && ictalMin < currMin
                currMin = ictalMin;
            end
            if dlrmBoolean == 1 && dlrmMin < currMin
                currMin = dlrmMin;
            end
            
            if ictalBoolean == 1 && ictalMax > currMax
                currMax = ictalMax;
            end
            if dlrmBoolean == 1 && dlrmMax > currMax
                currMax = dlrmMax;
            end
            xlim([currMin currMax]);
            set(gca, 'XScale', 'log');
            ylim([currMin currMax]);
            set(gca, 'YScale', 'log');


            axis square;
            plot([currMin currMax], [currMin currMax], 'b--', 'LineWidth', 1);
            xlabel("Control Spectral Power");
            ylabel("Ictal Spectral Power");
            title(currPatient);
            
            
            %Make a y = x line too...and separate the ictal and delirium
            %because you can't make shit out of it when they're all in one
            %tiny subplot!!!!!!!!!!!!!!!!
            
        end

        currPatient = newPatient;
        ictalBoolean = 0;
        dlrmBoolean = 0;
        patientCounter = patientCounter + 1;
    else 
        newPatientBoolean = 0;
    end
    
    if lastRunBoolean == 1
        plot(0, 0, 'Color', 'red');
        plot(0, 0, 'Color', 'green');
        legend("y = x Line", "Non-Delirious", "Delirious");
        sgtitle([band ' Band']);
        break;
    end
    % Make an rray of ROIs where each index matches the ROI's corresponding spectral power in
    % the x and y arrays
    currROIsCell = keys(dictionaryCells{iCombo, 2});
    %Have to go index by index to convert the cell of strings into an array
    currROIsArray = string(zeros(1, length(currROIsCell)));
    for i = 1:length(currROIsCell)
        currROIsArray(1, i) = string(currROIsCell{1, i});
    end


    % Make an x array of control spectral powers
    % Make two optional y arrays of ictal spectral powers
    currCondition = char(dictionaryCells{iCombo, 1});
    currCondition = currCondition(6:end);
    if strcmp(currCondition, "Ctrl") == 1
        ctrlSpectral = cell2mat(values(dictionaryCells{iCombo, 2}));
    elseif strcmp(currCondition, "Ictal") == 1
        ictalSpectral = cell2mat(values(dictionaryCells{iCombo, 2}));
        ictalBoolean = 1;
    elseif strcmp(currCondition, "IctalDlrm") == 1
        dlrmSpectral = cell2mat(values(dictionaryCells{iCombo, 2}));
        dlrmBoolean = 1;
    end

    
end

%Do another set of subplots of x = ictal non-delirious spectral power and y = ictal delirious
%spectral power
h = figure('Position',[65 158 1485 792]);
subplotsSize = 2;

usedPatientCounter = 0;
for iPatient = 1:length(patientList)
    patientRows = ib==iPatient;
    patientConditions = dataTable.mapCondition(patientRows);
    if length(patientConditions) ~= 3
        continue;
    else
        usedPatientCounter = usedPatientCounter + 1;
    end
    
    figure(h)
    subplot(subplotsSize,subplotsSize, usedPatientCounter);
    hold on;
    
    currDictCells = dictionaryCells(patientRows, 1:2);
    for iCondition = 1:length(currDictCells)
        currROIsCell = keys(currDictCells{iCondition, 2});
        %Have to go index by index to convert the cell of strings into an array
        currROIsArray = string(zeros(1, length(currROIsCell)));
        for i = 1:length(currROIsCell)
            currROIsArray(1, i) = string(currROIsCell{1, i});
        end

        currCondition = char(currDictCells{iCondition, 1});
        currCondition = currCondition(6:end);
        
        if strcmp(currCondition, "Ctrl") == 1
            ctrlSpectral = "Not Used";
        end
        
        if strcmp(currCondition, "Ictal") == 1
            ictalSpectral = cell2mat(values(currDictCells{iCondition, 2}));
        end

        if strcmp(currCondition, "IctalDlrm") == 1
            dlrmSpectral = cell2mat(values(currDictCells{iCondition, 2}));
        end
    end
    
    text(ictalSpectral, dlrmSpectral, currROIsArray, 'Color', 'blue', 'FontSize', 10, 'HorizontalAlignment','center');
  
    ictalMin = min(ictalSpectral);
    ictalMax = max(ictalSpectral);
    dlrmMin = min(dlrmSpectral);
    dlrmMax = max(dlrmSpectral);
    if ictalMin < dlrmMin
        currMin = ictalMin;
    else
        currMin = dlrmMin;
    end
    if ictalMax > dlrmMax
        currMax = ictalMax;
    else
        currMax = dlrmMax;
    end
    xlim([currMin currMax]);
    set(gca, 'XScale', 'log');
    ylim([currMin currMax]);
    set(gca, 'YScale', 'log');

    axis square;
    plot([currMin currMax], [currMin currMax], 'r--', 'LineWidth', 1);
    xlabel("Ictal Spectral Power");
    ylabel("IctalDlrm Spectral Power");
    title(patientList(iPatient));
    sgtitle([band ' Band']);
end

%% ROIs with largest changes in spectral power

topSpectralChange = string(zeros(usedPatientCounter, 11));

usedPatientCounter = 0;
for iPatient = 1:length(patientList)
    currPatient = patientList(iPatient);
    patientRows = ib==iPatient;

    patientConditions = dataTable.mapCondition(patientRows);
    if length(patientConditions) ~= 3
        continue;
    else
        usedPatientCounter = usedPatientCounter + 1;
    end
    
    
    currDictCells = dictionaryCells(patientRows, 1:2);
    ictalBoolean = 0;
    dlrmBoolean = 0;
    for iCondition = 1:length(currDictCells)
        currROIsCell = keys(currDictCells{iCondition, 2});
        %Have to go index by index to convert the cell of strings into an array
        currROIsArray = string(zeros(1, length(currROIsCell)));
        for i = 1:length(currROIsCell)
            currROIsArray(1, i) = string(currROIsCell{1, i});
        end

        currCondition = char(currDictCells{iCondition, 1});
        currCondition = currCondition(6:end);
        
        if strcmp(currCondition, "Ctrl") == 1
            ctrlSpectral = "Not Used";
        end
        if strcmp(currCondition, "Ictal") == 1
            ictalSpectral = cell2mat(values(currDictCells{iCondition, 2}));
            ictalBoolean = 1;
        end
        if strcmp(currCondition, "IctalDlrm") == 1
            dlrmSpectral = cell2mat(values(currDictCells{iCondition, 2}));
            dlrmBoolean = 1;
        end
            
        %After the last condition:
        if iCondition == length(currDictCells)
            if ictalBoolean == 1 && dlrmBoolean == 1
                % Get the change in spectral power and create a map from
                % the change to the corresponding ROI so that once you sort
                % the changes you can build back the top ROIs
                changeSpectral = abs(dlrmSpectral - ictalSpectral);
                change2roi = containers.Map(changeSpectral, currROIsArray);
                
                changeSorted = sort(changeSpectral, 'descend');

                topSpectralChange(usedPatientCounter, 1) = currPatient;
                for iVal = 1:10
                    currROI = string(change2roi(changeSorted(iVal)));
                    topSpectralChange(usedPatientCounter, iVal+1) = currROI;
                end
            end
        end

    end
    %End of Conditions Loop

end
%End of Patients Loop


% ROIs with top 10 large change per patient and put it into a cell
% --> First column patient name and the rest the ROIs

% Could do a dictionary with ROI key and value of number of times in the
% top 10 but there's only 3 subjects so not worth it yet...

%% Top Spectral Power Change Output

fileOutPath = ['C:\Users\Benjamin Hayum\Desktop\spectralPowerROIList\DataTables\topSpectralChange_ICTALvsICTALDLRM\' band 'SpectralChange_ICTALvsICTALDLRM.mat'];
save(fileOutPath, 'topSpectralChange');

%whole function end    
end