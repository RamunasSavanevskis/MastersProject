classdef GroupedData
    properties  
        originalAlpha
        groupIndex
        data
        groupedAlpha
        strainRanges
        yPositions

        groupAlphaChanges
        disagreementData
        

    end
    methods
        function obj = GroupedData(groupIndex, data, strainRanges, groupedAlpha, yPositions)
            obj.groupIndex = groupIndex;
            obj.data = data;
            obj.strainRanges = flip(abs(strainRanges .')) .';
            obj.groupedAlpha = groupedAlpha;
            obj.originalAlpha = groupedAlpha(groupIndex);
            obj.yPositions = yPositions;
            obj.groupAlphaChanges = 0;
            obj.disagreementData = struct();
        end
        function stresses = noCCalcStress(obj, exx)
            [~, ~, binIndices] = histcounts(abs(exx), 'BinEdges', obj.strainRanges);  %find out which "group" each strain belongs to

            stresses = exx .* obj.groupedAlpha(binIndices);  %multiply each strain by appropriate constant
        end
        function obj = recursiveNoC(obj, precision, iteration)  %basic sxx = alpha * exx, recursing on itself to get increasingly close to a good alpha for the group
            if iteration > 1
                previousDisagreementData = obj.disagreementData(iteration-1);
                previousAverageDisagreement = previousDisagreementData.averageDisagreement;  %check previous to see how want to change
                previousChange = obj.groupAlphaChanges(iteration-1);
                if iteration > 2
                    prePreviousDisagreementData = obj.disagreementData(iteration-2);
                    prePreviousAverageDisagreement = prePreviousDisagreementData.averageDisagreement;
                    if previousAverageDisagreement / prePreviousAverageDisagreement < 0  
                        %when change sign, must be value in between - can
                        %narrow down range
                    
                        precision = precision / 2;
                    end
                end
                if previousAverageDisagreement > 0
                    %decrease if too much moment/stress
                    newChange = previousChange - precision;
                else
                    %increase if too little moment/stress
                    newChange = previousChange + precision;
                end
                obj.groupAlphaChanges(iteration) = newChange;
                obj.groupedAlpha(obj.groupIndex:end) = obj.originalAlpha + newChange;  %change the alpha constant 
            end
            
            disagreementVectorForIteration = zeros(length(obj.data),1);
            for i=1:length(obj.data)  %goes through all exx arrays in the group
                dataExamined = obj.data(i);
                exxTrying = dataExamined.exxDist;
                sxxTrying = obj.noCCalcStress(exxTrying);
                calcMoment = trapz(obj.yPositions, obj.yPositions.*sxxTrying) * 5 *10^-3;  %moment based on standard equation
                %calcForce = trapz(obj.yPositions,sxxTrying);   %use later
                basicMoment = dataExamined.load * dataExamined.position;  %just force * distance
                basicDisagreement = (calcMoment - basicMoment) / basicMoment;   %how much disagree check
                disagreementVectorForIteration(i,1) = basicDisagreement;
            end
            avgDisagreement = mean(disagreementVectorForIteration);
            obj.disagreementData(iteration,1).disagreementAll = disagreementVectorForIteration;  %keep track of disagreement distribution - relevant later
            obj.disagreementData(iteration,1).averageDisagreement = avgDisagreement;  %average is most basic check - if low average, then agree well so good stress
            if abs(avgDisagreement) < 5e-4   %once agreement is good enough, can finish
                return
            end
            iteration = iteration + 1;  
            obj = obj.recursiveNoC(precision,iteration);  %recurses back, trying again with different alpha constant
        end
    end
end

