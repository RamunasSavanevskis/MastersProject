classdef DataExtractor2
    properties
        loadsUsed
        directory     
        xPixelMap    
        yPixelMap    
        stepsArray    
        prelimData      
        yDispData
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
        %first and last image manually defined from convergence test
        startIndex = 97;  
        finishIndex = 681;  
    end
    methods
        function obj = DataExtractor2(directory)
            obj.directory = directory;
            
         
            obj.stepsArray = [obj.finishIndex - obj.startIndex + 1, 1];  %relevant files only examine                     
            obj = obj.extractAllData();  %extract raw data in basic form
            obj = obj.calibrate();  %image number -> loads, pixels -> positions
            obj = obj.allDataGen();  
            obj = obj.filterUnwanted(6,60,50);  %filter unclean data
            obj = obj.unravelAllData();
            
        end        
        function obj = extractAllData(obj)
            filePattern = fullfile(obj.directory, 'Ti_bend1_*.mat');  
            fileList = dir(filePattern);
            
            firstFileName = fileList(1).name;
            firstOut = load(fullfile(obj.directory, firstFileName));   %x and y values same for all - can use to map using this. also initialising
            obj.xPixelMap = firstOut.x(1, 1:end).';  
            obj.yPixelMap = firstOut.y(1:end, 1);  
            obj.prelimData = cell(obj.finishIndex-obj.startIndex,length(obj.xPixelMap)); 
            obj.yDispData = cell(obj.finishIndex-obj.startIndex,length(obj.xPixelMap));
            for i = obj.startIndex:obj.finishIndex   %loops through all relevant files - some images with no load, some post-load
                fileName = fileList(i).name;
                loadStr = regexp(fileName, '\d+', 'match');
                number = str2double(loadStr{end});
                obj.stepsArray(i - obj.startIndex + 1, 1) = number; 
                out1 = load(fullfile(obj.directory, fileName));
                exx = out1.exx;  
                v = out1.v;
                for xIndex = 1:length(obj.xPixelMap)
                    obj.prelimData{i - obj.startIndex + 1, xIndex} = exx(:, xIndex);  %strain data
                    obj.yDispData{i - obj.startIndex + 1, xIndex} = v(:, xIndex);  %y displacement - for aligning
                    %row = load, column = x position. Raw, unfiltered data

                end
            end
        end

        function obj = calibrate(obj)   %gets load and positions array from file number and x and y pixel
            obj.xPositions = obj.xPixelMap * ((23.5*10^-3)/(1900)) + (1.5*10^-3);
            obj.yPositions = ((41/6) - (obj.yPixelMap/63)) * 10^-3;
            obj.loads = zeros(length(obj.stepsArray),1);
            loadStep = 0;
            for i=1:length(obj.stepsArray)
                obj.loads(i,1) = loadStep;
                loadStep = loadStep + 2.5;  %going up 2.5N per image
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
                    obj.allData(r,c).exxDist = exxArr; %most important, raw part
                    %min and max for asymmetry analysis
                    obj.allData(r,c).maxExx = max(exxArr);  
                    obj.allData(r,c).minExx = min(exxArr);  
                    obj.allData(r,c).maxAbsExx = max(abs(exxArr));
                    %calibrated data put into allData
                    obj.allData(r,c).load = obj.loads(r);  
                    obj.allData(r,c).position = obj.xPositions(c);
                    obj.allData(r,c).moment = obj.xPositions(c) .* obj.loads(r);

                end
            end
            
        end
        function obj = filterUnwanted(obj, firstX, lastX, firstLoad)   %dont care about edges + low loads (just noise)
            %filter all 3 - data, x positions, loads
            obj.allDataFiltered = obj.allData(firstLoad:end,firstX:lastX);
            obj.xPositionsFiltered = obj.xPositions(firstX:lastX);
            obj.loadsFiltered = obj.loads(firstLoad:end);
        end
        function obj = unravelAllData(obj)  %to change format/type of data for easier future use
            obj.allDataUnordered = reshape(obj.allDataFiltered, [], 1);  %"unravelling" allData into Ax1 array           
            dataTable = struct2table(obj.allDataUnordered);
            sortedTable = sortrows(dataTable, 'moment', 'ascend');
            obj.allDataOrdered = table2struct(sortedTable);  %put in table format, in order of moment
        end
        function obj = groupBySpacing(obj, spacing)     %NO LONGER RELEVANT - DONE EXTERNALLY
            assignin("base", "sortedData", obj.allDataOrdered)
            maxValue = 0;  % going from 0
            minValue = -1*max([obj.allDataOrdered.maxAbsExx]);  % to most negative strain
            minExxValues = -1 * [obj.allDataOrdered.maxAbsExx];
            binEdges = minValue:spacing:maxValue;
            binEdges = [-inf binEdges 0];
            obj.strainRanges = binEdges .';
            %used to divide the data into chunks
            [~, ~, binIndices] = histcounts(minExxValues, 'BinEdges', binEdges);
            uniqueBins = unique(binIndices);
            %an array of structures with the same data but now split based on the lowest exx value
            
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
            axis equal;
            hold on;

            %plot the y = -x line for demonstration
            fplot(@(x) -x, [-0.04, 0], 'r--', 'LineWidth', 2);

            %for demonstration purposes of zooming in
            fplot(@(x) x+0.025, [-0.04, 0], 'b--', 'LineWidth', 2);
            fplot(@(x) x+0.026, [-0.04, 0], 'b--', 'LineWidth', 2);

            title("pos vs neg - all scatter");
            set(gca, 'XDir', 'reverse');
            hold off;

            line1 = @(x) x + 0.025;   % First line (y = x + 0.02)
            line2 = @(x) x + 0.026;  % Second line (y = x + 0.021)

            %initialize a variable to hold the data points between the two lines
            pointsBetweenLines = [];

            %loop over each data point to find points between the two lines
            for i = 1:length(maxexx)
                x = minexx(i);    % Get x-coordinate from the negative strain data
                y = maxexx(i);    % Get y-coordinate from the positive strain data

                % Evaluate y-values of the two lines at the current x
                y1 = line1(x);
                y2 = line2(x);

                % Check if the y-value of the current data point is between the two lines
                if y > y1 && y < y2
                    pointsBetweenLines = [pointsBetweenLines; x, y]; % Add (x, y) to the result
                end
            end

            %display and plot points - again for demonstration
            figure;
            scatter(pointsBetweenLines(:, 1), pointsBetweenLines(:, 2));
            xMin = min(pointsBetweenLines(:, 1));
            xMax = max(pointsBetweenLines(:, 1));
            yMin = min(pointsBetweenLines(:, 2));
            yMax = max(pointsBetweenLines(:, 2));
            xlim([-1*max(abs(xMin),abs(yMin)),-1*min(abs(xMin),abs(yMin))])
            ylim([1*min(abs(xMin),abs(yMin)),1*max(abs(xMin),abs(yMin))])
            set(gca, 'XDir', 'reverse');
            axis equal
            hold on;
            fplot(@(x) -x, [min(pointsBetweenLines(:, 1)), max(pointsBetweenLines(:, 1))], 'r--', 'LineWidth', 2);
            axis equal
            title('Points Between Lines with y=-x Line');
            xlim([-1*max(abs(xMin),abs(yMin)),-1*min(abs(xMin),abs(yMin))])
            ylim([1*min(abs(xMin),abs(yMin)),1*max(abs(xMin),abs(yMin))])
            hold off;
        end

        
        function obj = individualExxGraph(obj, loadIndex, xPosIndex)
            data = obj.prelimData{loadIndex, xPosIndex};
            figure;
            yPositionVector = obj.yPositions + 2.5e-3;
            scatter(yPositionVector, data)
            hold on;
            p = polyfit(yPositionVector, data, 1);
            yFit = polyval(p, yPositionVector);
            plot(yPositionVector, yFit, '--', 'LineWidth', 1.5, 'Color', 'r');
            xAtY0 = -p(2)/p(1); % Solve p(1)*x + p(2) = 0
            disp(xAtY0)

            % Compute R^2 value
            yMean = mean(data);
            SS_tot = sum((data - yMean).^2);
            SS_res = sum((data - yFit).^2);
            R2 = 1 - (SS_res / SS_tot);

            % Add text to the plot
            xLimits = xlim;
            yLimits = ylim;

            % Add text to the top-right corner
            text(xLimits(2), yLimits(2), ...
                sprintf('y = %.4fx + %.4f\nR^2 = %.4f\nneutral axis position: %.4f', p(1), p(2), R2, xAtY0*1E3), ...
                'FontSize', 10, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'BackgroundColor', 'w');
            hold off
        end
    end
end
