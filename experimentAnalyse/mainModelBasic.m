clear
close all
clc

%obj = DataExtractor2("C:\Users\savan\Downloads\LimitedAOI6\LimitedAOI6");  %gathering all the data from the .mat files
%save('mySavedObject.mat', 'obj');
load('mySavedObject.mat', 'obj');
alpha1s = [];
alpha2s = [];
c1s = [];
c2s = [];
prevNegLimits = [];
prevPosLimits = [];

chunkSize = 1200;
numChunks = ceil(length(obj.allDataOrdered) / chunkSize);

groupSolvers = cell(numChunks, 1);

for i = 1:numChunks
    startIdx = (i-1) * chunkSize + 1;
    endIdx = min(i * chunkSize, length(obj.allDataOrdered));
    
    dataChunk = obj.allDataOrdered(startIdx:endIdx);
    
    groupSolvers{i} = groupSolver(dataChunk);
end
posRange = -100000e6:100000e6:100000e6;  %guesses
negRange = -100000e6:100000e6:100000e6;  
testRange = -20000E6:20000E6:20000E6;
testRange = testRange + [1000E6, 1000E6, 1000E6];
%{
firstOne = groupSolver(obj.allDataOrdered(1000:16000));
firstOne = firstOne.matrixParametersGet(1, alpha1s, alpha2s, c1s, c2s, prevNegLimits , prevPosLimits);

[correctPosStart, correctNegStart] = findBest(firstOne, 100E6,-100E6);
[firstOne, outAlpha1, outAlpha2, outC1, outC2, p1val, p2val] = firstOne.zoomInFixedIntercepts(correctPosStart, correctNegStart, testRange, testRange, 1,1,1);
fprintf('Highest positive: %.5g\n', firstOne.highestMaxExx);
fprintf('Highest negative: %.5g\n', firstOne.highestMinExx);
fprintf('Lowest positive: %.5g\n', firstOne.lowestMaxExx);
fprintf('Lowest negative: %.5g\n', firstOne.lowestMinExx);
fprintf('Pos start: %.5g\n', correctPosStart);
fprintf('Neg start %.5g\n', correctNegStart);
%}
%[firstOne, outAlpha1, outAlpha2, outC1, outC2, p1val, p2val, averageM, averageN] = firstOne.manualCheck(0,0, 112E9 * firstOne.highestMaxExx, 112E9 * firstOne.highestMinExx, 1);

%{
if 1 == 1
    alpha1s = [1.190579029063038e+11];
    alpha2s = [1.067265274467504e+11];
    c1s = [-1.810840271411017e+06];
    c2s = [9.654067665404540e+05];
else
    alpha1s = [1.13e11];
    alpha2s = [1.13e11];
    c1s = [0];
    c2s = [0];
end

prevNegLimits = [0];
prevPosLimits = [0];
firstOne = groupSolver(obj.allDataOrdered(4001:8000));
firstOne = firstOne.matrixParametersGet(0, alpha1s, alpha2s, c1s, c2s, prevNegLimits , prevPosLimits);

[correctPosStart, correctNegStart] = findBest(firstOne, 100E6,-100E6);
%[firstOne, outAlpha1, outAlpha2, outC1, outC2, p1val, p2val] = firstOne.zoomInFixedIntercepts(correctPosStart, correctNegStart, testRange, testRange, 1,1,1);
fprintf('Highest positive: %.5g\n', firstOne.highestMaxExx);
fprintf('Highest negative: %.5g\n', firstOne.highestMinExx);
fprintf('Lowest positive: %.5g\n', firstOne.lowestMaxExx);
fprintf('Lowest negative: %.5g\n', firstOne.lowestMinExx);

fprintf('Pos start: %.5g\n', correctPosStart);
fprintf('Neg start %.5g\n', correctNegStart);

%}
%{
alpha1s = [1.127350283239172e+11;1.190579029063038e+11];
alpha2s = [1.128468193643729e+11;1.067265274467504e+11];
c1s = [7.185940716737241e+06;-1.810840271411017e+06];
c2s = [1.288853934416172e+07;9.654067665404540e+05];

prevNegLimits = [-0.002092301902307;0];
prevPosLimits = [0.002006944651221;0];
firstOne = groupSolver(obj.allDataOrdered(8001:12000));
firstOne = firstOne.matrixParametersGet(0, alpha1s, alpha2s, c1s, c2s, prevNegLimits , prevPosLimits);

[correctPosStart, correctNegStart] = findBest(firstOne, 100E6,-100E6);
correctPosStart2 = mean(correctPosStart);
correctNegStart2 = mean(correctNegStart);
[firstOne, outAlpha1, outAlpha2, outC1, outC2, p1val, p2val] = firstOne.zoomInFixedIntercepts(correctPosStart2, correctNegStart2, testRange, testRange, 1,1,1);
fprintf('Highest positive: %.5g\n', firstOne.highestMaxExx);
fprintf('Highest negative: %.5g\n', firstOne.highestMinExx);
fprintf('Lowest positive: %.5g\n', firstOne.lowestMaxExx);
fprintf('Lowest negative: %.5g\n', firstOne.lowestMinExx);

fprintf('Pos start: %.5g\n', correctPosStart2);
fprintf('Neg start %.5g\n', correctNegStart2);
%}

alpha1s = [1.179454256732126e+11;1.127350283239172e+11;1.190579029063038e+11];
alpha2s = [9.309492245167656e+10;1.128468193643729e+11;1.067265274467504e+11];
c1s = [-2.535245506731653e+07;7.185940716737241e+06;-1.810840271411017e+06];
c2s = [-6.699774417417467e+07;1.288853934416172e+07;9.654067665404540e+05];

prevNegLimits = [-0.003335041479044;-0.002092301902307;0];
prevPosLimits = [0.003220048722116;0.002006944651221;0];
firstOne = groupSolver(obj.allDataOrdered(12001:16000));
firstOne = firstOne.matrixParametersGet(0, alpha1s, alpha2s, c1s, c2s, prevNegLimits , prevPosLimits);

[correctPosStart, correctNegStart] = findBest(firstOne, 1000E6,-200E6);
correctPosStart2 = mean(correctPosStart);
correctNegStart2 = mean(correctNegStart);
[firstOne, outAlpha1, outAlpha2, outC1, outC2, p1val, p2val] = firstOne.zoomInFixedIntercepts(correctPosStart2, correctNegStart2, testRange, testRange, 1,1,1);
fprintf('Highest positive: %.5g\n', firstOne.highestMaxExx);
fprintf('Highest negative: %.5g\n', firstOne.highestMinExx);
fprintf('Lowest positive: %.5g\n', firstOne.lowestMaxExx);
fprintf('Lowest negative: %.5g\n', firstOne.lowestMinExx);
alpha1s = [1.179454256732126e+11;1.127350283239172e+11;1.190579029063038e+11];
alpha2s = [9.309492245167656e+10;1.128468193643729e+11;1.067265274467504e+11];
c1s = [-2.535245506731653e+07;7.185940716737241e+06;-1.810840271411017e+06];
c2s = [-6.699774417417467e+07;1.288853934416172e+07;9.654067665404540e+05];

prevNegLimits = [-0.003335041479044;-0.002092301902307;0];
prevPosLimits = [0.003220048722116;0.002006944651221;0];
firstOne = groupSolver(obj.allDataOrdered(13001:16000));
firstOne = firstOne.matrixParametersGet(0, alpha1s, alpha2s, c1s, c2s, prevNegLimits , prevPosLimits);

[correctPosStart, correctNegStart] = findBest(firstOne, 1000E6,-200E6);
correctPosStart2 = mean(correctPosStart);
correctNegStart2 = mean(correctNegStart);
[firstOne, outAlpha1, outAlpha2, outC1, outC2, p1val, p2val] = firstOne.zoomInFixedIntercepts(correctPosStart2, correctNegStart2, testRange, testRange, 1,1,1);
fprintf('Highest positive: %.5g\n', firstOne.highestMaxExx);
fprintf('Highest negative: %.5g\n', firstOne.highestMinExx);
fprintf('Lowest positive: %.5g\n', firstOne.lowestMaxExx);
fprintf('Lowest negative: %.5g\n', firstOne.lowestMinExx);

fprintf('Pos start: %.5g\n', correctPosStart2);
fprintf('Neg start %.5g\n', correctNegStart2);
%}
%{
alpha1s = [9.052517358388567e+10;1.179454256732126e+11;1.127350283239172e+11;1.190579029063038e+11];
alpha2s = [1.083500982396138e+11;9.309492245167656e+10;1.128468193643729e+11;1.067265274467504e+11];
c1s = [9.446685493413323e+07;-2.535245506731653e+07;7.185940716737241e+06;-1.810840271411017e+06];
c2s = [-1.202507260847992e+07;-6.699774417417467e+07;1.288853934416172e+07;9.654067665404540e+05];

prevNegLimits = [-0.004569230773127;-0.003335041479044;-0.002092301902307;0];
prevPosLimits = [0.004468521908529;0.003220048722116;0.002006944651221;0];
firstOne = groupSolver(obj.allDataOrdered(16001:20000));
firstOne = firstOne.matrixParametersGet(0, alpha1s, alpha2s, c1s, c2s, prevNegLimits , prevPosLimits);

[correctPosStart, correctNegStart] = findBest(firstOne, 1000E6,-200E6);
correctPosStart2 = mean(correctPosStart);
correctNegStart2 = mean(correctNegStart);
[firstOne, outAlpha1, outAlpha2, outC1, outC2, p1val, p2val] = firstOne.zoomInFixedIntercepts(correctPosStart2, correctNegStart2, testRange, testRange, 1,1,1);
fprintf('Highest positive: %.5g\n', firstOne.highestMaxExx);
fprintf('Highest negative: %.5g\n', firstOne.highestMinExx);
fprintf('Lowest positive: %.5g\n', firstOne.lowestMaxExx);
fprintf('Lowest negative: %.5g\n', firstOne.lowestMinExx);

fprintf('Pos start: %.5g\n', correctPosStart2);
fprintf('Neg start %.5g\n', correctNegStart2);
%}
%{
alpha1s = [5.255346061793876e+10;9.052517358388567e+10;1.179454256732126e+11;1.127350283239172e+11;1.190579029063038e+11];
alpha2s = [7.832297845542671e+10;1.083500982396138e+11;9.309492245167656e+10;1.128468193643729e+11;1.067265274467504e+11];
c1s = [3.248770865414918e+08;9.446685493413323e+07;-2.535245506731653e+07;7.185940716737241e+06;-1.810840271411017e+06];
c2s = [-2.094866269280863e+08;-1.202507260847992e+07;-6.699774417417467e+07;1.288853934416172e+07;9.654067665404540e+05];

prevNegLimits = [-0.0058072;-0.004569230773127;-0.003335041479044;-0.002092301902307;0];
prevPosLimits = [0.0057467;0.004468521908529;0.003220048722116;0.002006944651221;0];
firstOne = groupSolver(obj.allDataOrdered(20001:24000));
firstOne = firstOne.matrixParametersGet(0, alpha1s, alpha2s, c1s, c2s, prevNegLimits , prevPosLimits);

[correctPosStart, correctNegStart] = findBest(firstOne, 1000E6,-200E6);
correctPosStart2 = mean(correctPosStart);
correctNegStart2 = mean(correctNegStart);
[firstOne, outAlpha1, outAlpha2, outC1, outC2, p1val, p2val] = firstOne.zoomInFixedIntercepts(correctPosStart2, correctNegStart2, testRange, testRange, 1,1,1);
fprintf('Highest positive: %.5g\n', firstOne.highestMaxExx);
fprintf('Highest negative: %.5g\n', firstOne.highestMinExx);
fprintf('Lowest positive: %.5g\n', firstOne.lowestMaxExx);
fprintf('Lowest negative: %.5g\n', firstOne.lowestMinExx);

fprintf('Pos start: %.5g\n', correctPosStart2);
fprintf('Neg start %.5g\n', correctNegStart2);
%}
%{
alpha1s = [3.365119163259831e+10;5.255346061793876e+10;9.052517358388567e+10;1.179454256732126e+11;1.127350283239172e+11;1.190579029063038e+11];
alpha2s = [-4.416371699320322e+10;7.832297845542671e+10;1.083500982396138e+11;9.309492245167656e+10;1.128468193643729e+11;1.067265274467504e+11];
c1s = [5.116214551637038e+08;3.248770865414918e+08;9.446685493413323e+07;-2.535245506731653e+07;7.185940716737241e+06;-1.810840271411017e+06];
c2s = [-1.190923116548040e+09;-2.094866269280863e+08;-1.202507260847992e+07;-6.699774417417467e+07;1.288853934416172e+07;9.654067665404540e+05];

prevNegLimits = [-0.0070798;-0.0058072;-0.004569230773127;-0.003335041479044;-0.002092301902307;0];
prevPosLimits = [0.0071562;0.0057467;0.004468521908529;0.003220048722116;0.002006944651221;0];
firstOne = groupSolver(obj.allDataOrdered(24001:end));
firstOne = firstOne.matrixParametersGet(0, alpha1s, alpha2s, c1s, c2s, prevNegLimits , prevPosLimits);

[correctPosStart, correctNegStart] = findBest(firstOne, 1000E6,-200E6);
correctPosStart2 = mean(correctPosStart);
correctNegStart2 = mean(correctNegStart);
[firstOne, outAlpha1, outAlpha2, outC1, outC2, p1val, p2val] = firstOne.zoomInFixedIntercepts(correctPosStart2, correctNegStart2, testRange, testRange, 1,1,1);
fprintf('Highest positive: %.5g\n', firstOne.highestMaxExx);
fprintf('Highest negative: %.5g\n', firstOne.highestMinExx);
fprintf('Lowest positive: %.5g\n', firstOne.lowestMaxExx);
fprintf('Lowest negative: %.5g\n', firstOne.lowestMinExx);

fprintf('Pos start: %.5g\n', correctPosStart2);
fprintf('Neg start %.5g\n', correctNegStart2);
%}
%{
alpha1s = [-2.616064117944069e+10;3.365119163259831e+10;5.255346061793876e+10;9.052517358388567e+10;1.179454256732126e+11;1.127350283239172e+11;1.190579029063038e+11];
alpha2s = [2.358355331018842e+10;-4.416371699320322e+10;7.832297845542671e+10;1.083500982396138e+11;9.309492245167656e+10;1.128468193643729e+11;1.067265274467504e+11];
c1s = [1.152110278990379e+09;5.116214551637038e+08;3.248770865414918e+08;9.446685493413323e+07;-2.535245506731653e+07;7.185940716737241e+06;-1.810840271411017e+06];
c2s = [-5.835324444823401e+08;-1.190923116548040e+09;-2.094866269280863e+08;-1.202507260847992e+07;-6.699774417417467e+07;1.288853934416172e+07;9.654067665404540e+05];

prevNegLimits = [-0.0084492;-0.0070798;-0.0058072;-0.004569230773127;-0.003335041479044;-0.002092301902307;0];
prevPosLimits = [0.0088703;0.0071562;0.0057467;0.004468521908529;0.003220048722116;0.002006944651221;0];
firstOne = groupSolver(obj.allDataOrdered(28001:end));
firstOne = firstOne.matrixParametersGet(0, alpha1s, alpha2s, c1s, c2s, prevNegLimits , prevPosLimits);

[correctPosStart, correctNegStart] = findBest(firstOne, 1000E6,-200E6);
correctPosStart2 = mean(correctPosStart);
correctNegStart2 = mean(correctNegStart);
[firstOne, outAlpha1, outAlpha2, outC1, outC2, p1val, p2val] = firstOne.zoomInFixedIntercepts(correctPosStart2, correctNegStart2, testRange, testRange, 1,1,1);
fprintf('Highest positive: %.5g\n', firstOne.highestMaxExx);
fprintf('Highest negative: %.5g\n', firstOne.highestMinExx);
fprintf('Lowest positive: %.5g\n', firstOne.lowestMaxExx);
fprintf('Lowest negative: %.5g\n', firstOne.lowestMinExx);

fprintf('Pos start: %.5g\n', correctPosStart2);
fprintf('Neg start %.5g\n', correctNegStart2);
%}
function [correctPosStart, correctNegStart] = findBest(group, posStartTry, negStartTry)
    posRange = -100000e6:100000e6:100000e6;  
    negRange = -100000e6:100000e6:100000e6;  %large range - covers as many possibilities as can
    
    %initial guess
    [group, ~,~,~,~, p1val, p2val] = group.zoomInFixedIntercepts(posStartTry, negStartTry, posRange, negRange, 1, 1, 0);
    primaryP1 = p1val;
    primaryP2 = p2val;
    
    %guesses make so can then find where 0
    posStarts = -400e6:20E6:1200E6;
    negStarts = -1200E6:20E6:400e6;
    p1_medians = zeros(length(posStarts), length(negStarts));
    p2_medians = zeros(length(posStarts), length(negStarts));
    for psi=1:length(posStarts)
        posStart = posStarts(psi);
        for nsi=1:length(negStarts)
            negStart = negStarts(nsi);

            if negStart == 0 && posStart == 0
                p1_medians(nsi,psi) = NaN;
                p2_medians(nsi,psi) = NaN;
                continue 
            end
            [group, ~,~,~,~, p1val, p2val] = group.zoomInFixedIntercepts(posStart, negStart, posRange, negRange, 1, 1, 0);
            p1Ratio = p1val ./ primaryP1;
            p2Ratio = p2val ./ primaryP2;
            p1_median = median(p1Ratio);
            p2_median = median(p2Ratio);
            p1_medians(nsi,psi) = p1_median;
            p2_medians(nsi,psi) = p2_median;

        end
    end
    [X,Y] = meshgrid(posStarts*1E-6, negStarts*1E-6);
    Z1 = p1_medians;
    Z2 = p2_medians;

    c1i = contourc(X(1,:), Y(:,1), Z1, [0 0]); 
    c2i = contourc(X(1,:), Y(:,1), Z2, [0 0]); %contour lines for where Zs are 0

    %extract (x, y) points for each contour
    [x1, y1] = extractContourPoints(c1i);
    [x2, y2] = extractContourPoints(c2i);

    %find intersections
    [xi, yi] = polyxpoly(x1, y1, x2, y2); % Find where the contours intersect
    correctPosStart = xi * 1e6;
    correctNegStart = yi * 1e6;
    %function to extract contour points
    function [x, y] = extractContourPoints(c)
        x = [];
        y = [];
        idx = 1;
        while idx < size(c, 2)
            n = c(2, idx); 
            x = [x, c(1, idx+1:idx+n)]; 
            y = [y, c(2, idx+1:idx+n)]; 
            idx = idx + n + 1;
        end
    end
    figure;
    surf(X,Y,Z1, 'EdgeColor', 'none');
    xlabel("Positive Start Stress (MPa)")
    ylabel("Negative Start Stress (MPa)")
    title("Moment Disagreement Distribution Scaling for Different Starting Stresses")
    figure;
    surf(X,Y,Z2, 'EdgeColor', 'none');
    xlabel("Positive Start Stress (MPa)")
    ylabel("Negative Start Stress (MPa)")
    title("Normal Distribution Scaling for Different Starting Stresses")
    
    figure; 
    hold on;
    plot(x1, y1, 'r', 'LineWidth', 2); %Z1 = 0 contour
    plot(x2, y2, 'b', 'LineWidth', 2); %Z2 = 0 contour
    plot(xi, yi, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k'); %intersections
    legend('Moment Disagreement Ratio = 0', 'Normal Ratio = 0');
    xlabel("Positive Start Stress (MPa)")
    ylabel("Negative Start Stress (MPa)")
    title("Contour Plot of Moment Disagreement and Normal Ratio Shrinking to 0")
    grid on;
    

end


