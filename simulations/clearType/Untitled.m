 filter = [2 6];
 %hold on;
 totAcc = zeros(5,14);
 totErr = zeros(5,14);
 for curIter = 1 : 10 
 for iFilter = -2 : 2
     kvDist = 1;
     for vDist = 0.7: 0.1 : 2
        defocus = iFilter;
        s_displayFilters;
        totAcc(iFilter+3, kvDist)  = totAcc(iFilter+3, kvDist) + acc;
        totErr(iFilter+3, kvDist) = totErr(iFilter+3, kvDist) + err;
        kvDist = kvDist + 1;
     end
     %errorbar(totAcc(iFilter-2,:), totErr(iFilter-2,:));
 end
 end