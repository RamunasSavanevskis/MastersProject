clear
close all
clc

obj = DataExtractor("C:\Users\mans4449\Desktop\DIC Tests part I\LimitedAOI6\");  %gathering all the data from the .mat files
obj = obj.groupBySpacing(0.001);  %group data in minimum strain intervals of 0.1%
groupedData = obj.groupedData;   %a structure array containing all the relevant data in groups based on mininimum exx
strainRanges = obj.strainRanges;
groupedAlpha = 120E9 * ones(length(strainRanges)+1, 1);  %generic starting point - different alpha for each
%with future ideas for implementation, may no longer be relevant - but
%keeping for now
yPositions = obj.yPositions;
for i=1:length(groupedData)
    newObj = GroupedData(groupedData(i).groupIndex,groupedData(i).data, strainRanges, groupedAlpha, yPositions);
    newObj = newObj.recursiveNoC(32E9,1);  %how far move the alpha (32GPa to start with). Iteration always = 1 at start
    display(i)  %just simple track of how long each takes
    groupedAlpha = newObj.groupedAlpha;  %with each group iteration, the groupedAlpha updates from the groupIndex to the end
    %groupIndex:end rather than just groupIndex purely for speed 
    groupedDataObjects(i,1) = newObj;  %stores all objects with alpha calculated in array for future analysis

end
final = newObj;  %using final because finished then
exxTry = linspace(-0.018,0.018, 1800) .';  %plot over some range - can go quite a bit further but interested here
sxxTry = final.noCCalcStress(exxTry);   %with the most recent alpha array, find what the stresses are
figure;
scatter(exxTry,sxxTry);  %graph without c constant


assignin("base", "allData", obj.allDataCell)