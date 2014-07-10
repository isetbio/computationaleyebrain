 filter = [2 3];
 %hold on;
 totAcc = zeros(4,14);
 totErr = zeros(4,14);
 for curIter = 1 : 10 
 for iFilter = 3 : 6
     kvDist = 1;
     for vDist = 0.7: 0.1 : 2
        humanConeDensities = [0 iFilter/10 0.9-iFilter/10 0.1];
        s_displayFilters_exp_coneDensity;
        totAcc(iFilter-2, kvDist)  = totAcc(iFilter-2, kvDist) + acc;
        totErr(iFilter-2, kvDist) = totErr(iFilter-2, kvDist) + err;
        kvDist = kvDist + 1;
     end
     %errorbar(totAcc(iFilter-2,:), totErr(iFilter-2,:));
 end
 end