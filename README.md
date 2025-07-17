Note: the code here was developed for my MEng project, and was developed just for my own use, meaning there are parts that need to be manually changed and are generally outdated. 
The basic idea of the useage is going into mainModelBasic and changing the commented-out blocks such as 

alpha1s = [1.127350283239172e+11;1.190579029063038e+11];
alpha2s = [1.128468193643729e+11;1.067265274467504e+11];
c1s = [7.185940716737241e+06;-1.810840271411017e+06];
c2s = [1.288853934416172e+07;9.654067665404540e+05];

prevNegLimits = [-0.002092301902307;0];
prevPosLimits = [0.002006944651221;0];
firstOne = groupSolver(obj.allDataOrdered(8001:12000));

to match the previous groups' alpha and c outputs
