classdef DataExtractor
    properties
        directory     
        xPixelMap    
        yPixelMap    
        stepsArray    
        prelimData      
        allDataCell
        xPositions
        yPositions
        loads
        allData

        allDataFiltered
        xPositionsFiltered
        loadsFiltered

        allDataUnordered
        allDataOrdered
        groupedData

        strainRanges
    end
    
    properties (Constant)
        %these 2 are right now manually defined - can add function to
        %examine the bend test csv data to align, but takes like 5 minutes
        %and procedure will be different once creep is there
        startIndex = 97;  
        finishIndex = 681;  
    end
    methods
        % Constructor to initialize the object and run the data extraction
        function obj = DataExtractor(directory)
            obj.directory = directory;
            
         
            obj.stepsArray = [obj.finishIndex - obj.startIndex + 1, 1];  %basically just file number, but corrected to be the relevant files                        
            obj = obj.extractAllData();
            obj = obj.calibrate();
            obj = obj.allDataGen();
            obj = obj.filterUnwanted(6,60,50);
            obj = obj.unravelAllData();
            
        end        
        % Method to extract data from .mat files
        function obj = extractAllData(obj)
            filePattern = fullfile(obj.directory, 'Ti_bend1_*.mat');   %change to input later - formatting of file
            fileList = dir(filePattern);
            
            firstFileName = fileList(1).name;
            firstOut = load(fullfile(obj.directory, firstFileName));   %x and y values same for all - can use to map using this. also initialising
            obj.xPixelMap = firstOut.x(1, 1:end).';  
            obj.yPixelMap = firstOut.y(1:end, 1);
            obj.prelimData = cell(obj.finishIndex-obj.startIndex,length(obj.xPixelMap));  
            for i = obj.startIndex:obj.finishIndex   %loops through all RELEVANT files - some images with no load, some post-load
                fileName = fileList(i).name;
                loadStr = regexp(fileName, '\d+', 'match');
                number = str2double(loadStr{end});
                obj.stepsArray(i - obj.startIndex + 1, 1) = number;  % Store the step number

                out1 = load(fullfile(obj.directory, fileName));
                exx = out1.exx;  

                
                for xIndex = 1:length(obj.xPixelMap)
                    obj.prelimData{i - obj.startIndex + 1, xIndex} = exx(:, xIndex);
                    %row = load, column = x position. Raw, unfiltered data
                    %each place in prelimData is a vector of the exx values
                    %across the y-axis for that x position + load
                end
            end
        end

        function obj = calibrate(obj)   %gets load and positions array from file number and x and y pixel
            obj.xPositions = obj.xPixelMap * ((23.5*10^-3)/(1900)) + (1.5*10^-3);
            obj.yPositions = ((41/6) - (obj.yPixelMap/63)) * 10^-3;
            %positions from pixels will be unnecessary when i use
            %calibrated data - this is very much "eye-balling" it
            obj.loads = zeros(length(obj.stepsArray),1);
            loadStep = 0;
            for i=1:length(obj.stepsArray)
                obj.loads(i,1) = loadStep;
                loadStep = loadStep + 2.5;
            end
        end
        function obj = allDataGen(obj)  %prelim data is just current relevant part (exx)
            [rows,columns] = size(obj.prelimData); 
            i = 0;
            obj.allData = struct();   %allData has same data as prelimData but is a structure with extra useful info
            for r=1:rows
                for c=1:columns
                    i = i + 1;
                    exxArr = obj.prelimData{r,c};
                    obj.allData(r,c).exxDist = exxArr;
                    obj.allData(r,c).maxExx = max(exxArr);  %might be used later for asymmetry analysis
                    obj.allData(r,c).minExx = min(exxArr);  %used for grouping - crucial
                    obj.allData(r,c).load = obj.loads(r);  
                    obj.allData(r,c).position = obj.xPositions(c);

                end
            end
            assignin("base", "alldata", obj.allData)
            
        end
        function obj = filterUnwanted(obj, firstX, lastX, firstLoad)   %dont care about edges + low loads (just noise)
            %firstX, lastX and firstLoad currently manual inputs from
            %looking at colour plots before - can maybe automate it, or
            %keep it manual
            [rows, cols] = size(obj.allData);
            disp(['Rows: ', num2str(rows), ', Columns: ', num2str(cols)]);

            obj.allDataFiltered = obj.allData(firstLoad:end,firstX:lastX);
            obj.xPositionsFiltered = obj.xPositions(firstX:lastX);
            obj.loadsFiltered = obj.loads(firstLoad:end);
        end
        function obj = unravelAllData(obj)
            obj.allDataUnordered = reshape(obj.allDataFiltered, [], 1);  %"unravelling" allData into Ax1 array 
            %allData technically unnecessary as AxB cell isn't necessary
            %for anything, but can be useful in some cases if you want to
            %explore whole field - either way works in graphing but for
            %manual inspection nice
                
            dataTable = struct2table(obj.allDataUnordered);
            sortedTable = sortrows(dataTable, 'minExx', 'descend');
            obj.allDataOrdered = table2struct(sortedTable);
        end
        function obj = groupBySpacing(obj, spacing)    
            assignin("base", "sortedData", obj.allDataOrdered)
            maxValue = 0;  % going from 0
            minValue = min([obj.allDataOrdered.minExx]);  % to most negative strain
            minExxValues = [obj.allDataOrdered.minExx];
            binEdges = minValue:spacing:maxValue;  
            binEdges = [-inf binEdges 0];
            obj.strainRanges = binEdges .';
            %used to divide the data into chunks
            [~, ~, binIndices] = histcounts(minExxValues, 'BinEdges', binEdges);  
            uniqueBins = unique(binIndices);
            %an array of structures with the same data but now split based
            %on the lowest exx value
            groupedStructs = struct('groupIndex', {}, 'data', {}); 
            
            for i = 1:length(uniqueBins)
                currentBin = uniqueBins(i);
                if currentBin == 0
                    continue; 
                end
                binData = obj.allDataOrdered(binIndices == currentBin);
                groupedStructs(end+1).groupIndex = max(uniqueBins) - currentBin + 1; 
                groupedStructs(end).data = binData; 
            end
        obj.groupedData = groupedStructs;
        obj.groupedData = obj.groupedData(end:-1:1);
        
        end

            
        function obj = minVMaxFig(obj)    %plot of the minimum and maximum exx for each load and position
            minexx = [obj.allDataOrdered.minExx] .';
            maxexx = [obj.allDataOrdered.maxExx] .';
            assignin("base", "minexx",minexx);
            assignin("base", "maxexx",maxexx);
            scatter(minexx, maxexx);
            set (gca, 'xdir', 'reverse' )
        end
        function obj = groupByAmount(obj,numGroups)   %implement later (not really necessary - effectively same previous

        end
    end
end
