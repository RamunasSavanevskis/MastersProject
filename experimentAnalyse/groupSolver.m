classdef groupSolver
    properties
        DataChunk
        maxAbsExx
        momentMatrix
        normalsMatrix
        momentVector




        prevNormalsVector
        prevMoments

        originalMomentVector

        lowestMaxExx
        lowestMinExx

        highestMaxExx 
        highestMinExx 
    end

    methods
        function obj = groupSolver(dataChunk)
            obj.DataChunk = dataChunk;
        end
        function obj = matrixParametersGet(obj, first, alpha1s, alpha2s, c1s, c2s, prevNegLimits , prevPosLimits)
            rawData = obj.DataChunk;
            % Extract values from the struct array
            momentValues = [rawData.moment];  % X-axis
            maxExxValues = [rawData.maxExx];  % Y-axis for first fit
            minExxValues = [rawData.minExx];  % Y-axis for second fit

            % fit cubic polynomial - generally good measure
            p1 = polyfit(momentValues, maxExxValues, 2);
            p2 = polyfit(momentValues, minExxValues, 2);  

            % get fitted results
            fittedMaxExx = polyval(p1, momentValues);
            fittedMinExx = polyval(p2, momentValues);

            obj.lowestMaxExx = fittedMaxExx(1);
            obj.lowestMinExx = fittedMinExx(1);
            obj.highestMaxExx = fittedMaxExx(end);
            obj.highestMinExx = fittedMinExx(end);
            if first == 1  %real data lowest wont be 0, but need to treat as is for maths
                obj.lowestMaxExx = 0;
                obj.lowestMinExx = 0;
            end
            obj.maxAbsExx = [rawData.maxAbsExx];
            yPos = linspace(0,5e-3,27);   %0 to 5mm, 27 data points along length. can change for different data
            obj.momentMatrix = zeros(length(rawData),4);
            obj.normalsMatrix = zeros(length(rawData),4);
            obj.momentVector = zeros(length(rawData),1);
            obj.originalMomentVector = zeros(length(rawData),1);

            huppers = zeros(length(rawData),1);
            obj.prevMoments = huppers;
            obj.prevNormalsVector = huppers;
            for i=1:length(rawData)
                
                moment = rawData(i).moment;
                obj.momentVector(i,1) = moment;
                obj.originalMomentVector(i,1) = moment;
                strainDistribution = flip(rawData(i).exxDist);
                
                distanceEquation = polyfit(yPos, strainDistribution, 1);
                NA = -distanceEquation(2) / distanceEquation(1);  %from y = mx + c
                hlower = 5e-3 - NA;
                hupper = -1 * NA;

                beta = distanceEquation(1);
                if first == 1   %works same way as if first != 1, but old, simple code
                   
                    af = (beta * hupper^3)/3;
                    bf = (hupper^2)/2;
                    ff = (beta * hlower^3)/3;
                    gf = (hlower^2)/2;
                    hf = (beta * hupper^2)/2;
                    kf = hupper - 0;
                    lf = (beta * hlower^2)/2;
                    mf = hlower - 0;
                    obj.momentMatrix(i,1) = af;
                    obj.momentMatrix(i,2) = bf;
                    obj.momentMatrix(i,3) = -ff;
                    obj.momentMatrix(i,4) = -gf;
                    obj.normalsMatrix(i,1) = hf;
                    obj.normalsMatrix(i,2) = kf;
                    obj.normalsMatrix(i,3) = -lf;
                    obj.normalsMatrix(i,4) = -mf;
                else
                    lowerLimits = [beta*hlower;obj.lowestMinExx;prevNegLimits];
                    upperLimits = [beta*hupper;obj.lowestMaxExx;prevPosLimits];                   
                    hupperPrevs = upperLimits / beta;
                    hlowerPrevs = lowerLimits / beta;
        
                    cubeLowerLimits = hlowerPrevs  .^3;
                    cubeUpperLimits = hupperPrevs .^3;
                    squareLowerLimits = hlowerPrevs  .^2;
                    squareUpperLimits = hupperPrevs .^2;

                    diffCubeUpper = cubeUpperLimits(1:end-1) - cubeUpperLimits(2:end);
                    diffSquareUpper = squareUpperLimits(1:end-1) - squareUpperLimits(2:end);
                    diffCubeLower = cubeLowerLimits(1:end-1) - cubeLowerLimits(2:end);
                    diffSquareLower = squareLowerLimits(1:end-1) - squareLowerLimits(2:end);
                    
                    A = zeros(length(diffCubeUpper), 4);
                    A(:, 1) = diffCubeUpper;  
                    A(:, 2) = diffSquareUpper;  
                    A(:, 3) = -diffCubeLower; 
                    A(:, 4) = -diffSquareLower;      
                    [nRows, nCols] = size(A);  %same dimensions
                    B = zeros(nRows, nCols);
                    B(:, 1:2:end) = beta / 3;  
                    B(:, 2:2:end) = 1 / 2;   
                    fullMomentMatrix = (A .* B);

                    diffSquareUpper = squareUpperLimits(1:end-1) - squareUpperLimits(2:end);
                    diffSquareLower = squareLowerLimits(1:end-1) - squareLowerLimits(2:end);
                    diffBaseUpper = hupperPrevs(1:end-1) - hupperPrevs(2:end);
                    diffBaseLower = hlowerPrevs(1:end-1) - hlowerPrevs(2:end);

                    A2 = zeros(nRows, nCols);  
                    A2(:, 1) = diffSquareUpper;
                    A2(:, 2) = diffBaseUpper;
                    A2(:, 3) = -diffSquareLower; 
                    A2(:, 4) = -diffBaseLower;  
                    [nRows, nCols] = size(A2);              
                    B2 = zeros(nRows, nCols);
                    B2(:, 1:2:end) = beta / 2; 
                    B2(:, 2:2:end) = 1;  
                    fullNormalMatrix = (A2 .* B2);
                    alphaCMatrix = zeros(4, length(alpha1s));
                    alphaCMatrix(1, 1:end) = alpha1s .';  
                    alphaCMatrix(2, 1:end) = c1s.';  
                    alphaCMatrix(3, 1:end) = alpha2s.';      
                    alphaCMatrix(4, 1:end) = c2s.';

                    obj.normalsMatrix(i,:) = fullNormalMatrix(1, :);
                    obj.momentMatrix(i,:) = fullMomentMatrix(1, :);
                    %sum diags for the maths to work
                    C = diag(fullMomentMatrix(2:end, :) * alphaCMatrix);
                    D = diag(fullNormalMatrix(2:end, :) * alphaCMatrix);
                    obj.prevNormalsVector(i) = sum(D);
                    obj.prevMoments(i,1) = sum(C);
                end
            end
            obj.momentMatrix = obj.momentMatrix * 5E-3;
            obj.normalsMatrix = obj.normalsMatrix * 5E-3;
            obj.prevMoments = obj.prevMoments * 5E-3;
            obj.prevNormalsVector = obj.prevNormalsVector * 5E-3;
            obj.momentVector = obj.momentVector - obj.prevMoments;
        end
        function [obj, outAlpha1, outAlpha2, outC1, outC2, p1val, p2val, averageM, averageN] = manualCheck(obj, posStart, negStart, posEnd, negEnd, show)
            startMinExx = obj.lowestMinExx;
            endMinExx = obj.highestMinExx;
            startMaxExx = obj.lowestMaxExx;
            endMaxExx = obj.highestMaxExx;
            outAlpha1 = (posEnd-posStart)/(endMaxExx-startMaxExx);
            outAlpha2 = (negEnd-negStart)/(endMinExx-startMinExx);  %CHANGE 0 LATER
            outC1 = posStart - (startMaxExx * outAlpha1);
            outC2 = negStart - (startMinExx * outAlpha2);

            resultMomentDifference = (obj.momentMatrix * [outAlpha1; outC1; outAlpha2; outC2] - obj.momentVector);
            resultNormals = (obj.normalsMatrix * [outAlpha1; outC1;outAlpha2; outC2]);

            p1 = polyfit(obj.momentVector, resultMomentDifference,3);
            p2 = polyfit(obj.momentVector, resultNormals,3);

            p1val = polyval(p1,obj.momentVector);
            p2val = polyval(p2,obj.momentVector);



            averageM = mean(resultMomentDifference);
            averageN = mean(resultNormals);

            if show == 1
                
                a1 = p1(1); b1 = p1(2); c1 = p1(3);
                a2 = p2(1); b2 = p2(2); c2 = p2(3);

                p1_deriv_eval = (a1 * 3) * (obj.momentVector .^ 2) + (b1 * 2) * obj.momentVector + c1;
                p2_deriv_eval = (a2 * 3) * (obj.momentVector .^ 2) + (b2 * 2) * obj.momentVector + c2;
                

                
                figure;
                scatter(obj.momentVector, resultMomentDifference);
                hold on;
                plot(obj.momentVector, p1val);
                title("Moment Disagreement Distribution")
                xlabel('Excess Moment (Nm)');
                ylabel('Moment disagreement (Nm)');
                ylim([-0.4; 0.4]);
                hold off;

                figure;
                scatter(obj.momentVector, resultNormals);
                hold on;
                plot(obj.momentVector, p2val);

                xlabel('Excess Moment (Nm)');
                ylabel('Normal (N)');
                title('Normal Distribution');
                hold off;

            end


        end
        function [obj, outAlpha1, outAlpha2, outC1, outC2, p1val, p2val] = zoomInFixedIntercepts(obj, posStart, negStart, posRange, negRange, posStartIndex, negStartIndex, show)
            %zoom-in algorithm to find where both go toward 0 - simple
            %planes for both
            endPosValues =  posRange;   
            endNegValues =  negRange;  
            startMinExx = obj.lowestMinExx;
            endMinExx = obj.highestMinExx;
            startMaxExx = obj.lowestMaxExx;
            endMaxExx = obj.highestMaxExx;
            if 1 == 2   %changing manually - could have input but no need
                startMinExx = 0;
                startMaxExx = 0;
            end
            stepPos = diff(endPosValues(1:2)) ;
            stepNeg = diff(endNegValues(1:2)) ;
            if abs(stepPos) <= 1E6 && abs(stepNeg) <= 1E6  %finish when both at least within 1MPa accuracy
               %find alpha,c values from limits
               posEnd = endPosValues(2);
               negEnd = endNegValues(2);
               outAlpha1 = (posEnd-posStart)/(endMaxExx-startMaxExx);
               outAlpha2 = (negEnd-negStart)/(endMinExx-startMinExx);
               outC1 = posStart - (startMaxExx * outAlpha1);
               outC2 = negStart - (startMinExx * outAlpha2);


               resultMomentDifference = (obj.momentMatrix * [outAlpha1; outC1; outAlpha2; outC2] - obj.momentVector);
               resultNormals = obj.normalsMatrix * [outAlpha1; outC1;outAlpha2; outC2] + obj.prevNormalsVector;

               
                p1 = polyfit(obj.momentVector, resultMomentDifference,3);
                p2 = polyfit(obj.momentVector, resultNormals,3);
                p1val = polyval(p1,obj.momentVector);
                p2val = polyval(p2,obj.momentVector);


               if show == 1  %show graphs when wanted

                   posLimits = linspace(startMaxExx, endMaxExx, 100);
                   negLimits = linspace(startMinExx, endMinExx, 100);
                   posStresses = (posLimits * outAlpha1) + outC1;
                   negStresses = (negLimits * outAlpha2) + outC2;
                   figure;
                   xlabel("Strain")
                   ylabel("Stress (MPa)")
                   title("Second Stress-Strain Section Found - Using Symmetric Previous Input")
                   hold on
                   h1 = plot(posLimits, 1E-6 * posStresses, 'b');
                   h2 = plot(negLimits, 1E-6 * negStresses, 'r');

                   %extract values (convert stress to MPa)
                   psTrainStart = posLimits(1);
                   psTrainEnd   = posLimits(end);
                   psStressStart = posStresses(1) * 1e-6;
                   psStressEnd   = posStresses(end) * 1e-6;

                   ngTrainStart = negLimits(1);
                   ngTrainEnd   = negLimits(end);
                   ngStressStart = negStresses(1) * 1e-6;
                   ngStressEnd   = negStresses(end) * 1e-6;

                   %plot invisible dummy points for legend handles
                   h3 = plot(NaN, NaN, 'b'); % Dummy positive start
                   h4 = plot(NaN, NaN, 'b'); % Dummy positive end
                   h5 = plot(NaN, NaN, 'r'); % Dummy negative start
                   h6 = plot(NaN, NaN, 'r'); % Dummy negative end

                   %create legend strings with strain and stress
                   legendLines = {
                       sprintf('Positive start: (%.4f, %.2f MPa)', psTrainStart, psStressStart)
                       sprintf('Positive end:   (%.4f, %.2f MPa)', psTrainEnd,   psStressEnd)
                       sprintf('Negative start: (%.4f, %.2f MPa)', ngTrainStart, ngStressStart)
                       sprintf('Negative end:   (%.4f, %.2f MPa)', ngTrainEnd,   ngStressEnd)
                       };
                   legend([h3 h4 h5 h6], legendLines, 'Location', 'northeast')
                   hold off
                   %extract coefficients from polyfit results
                   a1 = p1(1); b1 = p1(2); c1 = p1(3);
                   a2 = p2(1); b2 = p2(2); c2 = p2(3);
                   %compute derivatives manually
                   p1_deriv_eval = (a1 * 3) * (obj.momentVector .^ 2) + (b1 * 2) * obj.momentVector + c1;
                   p2_deriv_eval = (a2 * 3) * (obj.momentVector .^ 2) + (b2 * 2) * obj.momentVector + c2;  %both no longer needed but keeping
                   figure;
                   scatter(obj.momentVector, resultMomentDifference);
                   hold on;
                   plot(obj.momentVector, p1val);
                   title("Moment Disagreement Distribution")
                   xlabel('Excess Moment (Nm)');
                   ylabel('Moment disagreement (Nm)');

                   hold off;

                   figure;
                   scatter(obj.momentVector, resultNormals);
                   hold on;
                   plot(obj.momentVector, p2val);

                   xlabel('Excess Moment (Nm)');
                   ylabel('Normal (N)');
                   title('Normal Distribution');
                   hold off;                 
               end                
               return
            end
            numPosEnd = length(endPosValues);
            numNegEnd = length(endNegValues);
            averageMomentDistagreements = zeros(numNegEnd, numPosEnd);
            averageNormals = zeros(numNegEnd, numPosEnd);
            for i = 1:numPosEnd
                posEnd = endPosValues(i);
                for j = 1:numNegEnd
                    negEnd = endNegValues(j);
                    alpha1 = (posEnd-posStart)/(endMaxExx-startMaxExx);
                    alpha2 = (negEnd-negStart)/(endMinExx-startMinExx);
                    c1 = posStart - (startMaxExx * alpha1);
                    c2 = negStart - (startMinExx * alpha2);
                    resultMomentDifference = (obj.momentMatrix * [alpha1; c1; alpha2; c2] - obj.momentVector);
                    resultNormals = obj.normalsMatrix * [alpha1; c1; alpha2; c2] + obj.prevNormalsVector;
                    averageMomentDistagreements(i,j) = mean(resultMomentDifference);
                    averageNormals(i,j) = mean(resultNormals);
                end
            end
            
            Z1 = averageMomentDistagreements;
            Z2 = averageNormals;
            %{
            [X, Y] = meshgrid(endPosValues * 1E-6, endNegValues * 1E-6); 
            figure;
            surf(X,Y,Z1, 'EdgeColor','none')
            view(0,90)
            title("moment")
            xlabel("pos stress")
            ylabel("neg stress")
            figure;
            surf(X,Y,Z2, 'EdgeColor','none')
            view(0,90)
            title("normal")
            xlabel("pos stress")
            ylabel("neg stress")
            %}
            %here use 2x2 matrices representing corners of 3x3 signs (2x2
            %showing sign change or no) to find overlapping sign changes
            %sign change for both means both contain 0
            Z1Signs = flipud(sign(Z1));
            Z2Signs = flipud(sign(Z2));
            bothSigns = {Z1Signs, Z2Signs};
            bothValidCorners = cell(2,1);
            for i = 1:2
                signs = bothSigns{i};
                TL = signs(1:2,1:2);
                TR = signs(1:2,2:3);
                BL = signs(2:3,1:2);
                BR = signs(2:3,2:3);
                cornerValues = ones(2,2);
                if abs(sum(TL, 'all')) == 4
                    cornerValues(1,1) = 0;
                end

                if abs(sum(TR, 'all')) == 4
                    cornerValues(1,2) = 0;
                end

                if abs(sum(BL, 'all')) == 4
                    cornerValues(2,1) = 0;
                end

                if abs(sum(BR, 'all')) == 4
                    cornerValues(2,2) = 0;
                end
                bothValidCorners{i} = cornerValues;

            end
            overlap = bothValidCorners{1} .* bothValidCorners{2};
            if overlap == [1 0; 0 0]
                validNegRange = linspace(negRange(1), negRange(2),3);
                validPosRange = linspace(posRange(2), posRange(3), 3);
            elseif overlap == [0 1; 0 0]
                validNegRange = linspace(negRange(2), negRange(3),3);
                validPosRange = linspace(posRange(2), posRange(3), 3);
            elseif overlap == [0 0; 1 0]
                validNegRange = linspace(negRange(1), negRange(2),3);
                validPosRange = linspace(posRange(1), posRange(2), 3);
            elseif overlap == [0 0; 0 1]
                validNegRange = linspace(negRange(2), negRange(3),3);
                validPosRange = linspace(posRange(1), posRange(2), 3);
            elseif overlap == [1 1; 0 0]
                validNegRange = negRange;
                validPosRange = linspace(posRange(2), posRange(3), 3);
            elseif overlap == [1 0; 1 0]
                validNegRange = linspace(negRange(1), negRange(2),3);
                validPosRange = posRange;
            elseif overlap == [0 0; 1 1]
                validNegRange = negRange;
                validPosRange = linspace(posRange(1), posRange(2), 3);
            elseif overlap == [0 1; 0 1]
                validNegRange = linspace(negRange(2), negRange(3),3);
                validPosRange = posRange;
            end
            p1val = [];
            p2val = [];
            %iterate to "zoom-in" towards where both have sign change
            %(approaching 0)
            [obj, outAlpha1, outAlpha2, outC1, outC2, p1val, p2val] = obj.zoomInFixedIntercepts(posStart, negStart, validPosRange, validNegRange,  posStartIndex, negStartIndex, show);

           
        end
   
  
   
    end
end
