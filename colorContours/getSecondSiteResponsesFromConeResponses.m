function [blankSecondSiteResponses,testSecondSiteResponses] = getSecondSiteResponsesFromConeResponses(blankConeResponses,testConeResponses,theParams)
% [blankSecondSiteResponses,testSecondSiteResponses] = getSecondSiteResponsesFromConeResponses(blankConeResponses,testConeResponses,theParams)
%
% Get a set of second site responses from the cone responses.  How the second
% site responses are computed depends both on the type of the coneResponses
% structure, and on variables specified in the fields of theParams.
%
% We do both blank and test in same routine, so that the spatial draw of cones is the 
% same for each.  This is important.  Using different cones after the Poisson noise
% at the first stage can impose enough structure on the output that an svm can
% tell test from blank even with zero signal!  [We learned this the hard  way.]
%
% In the longer run, this needs to take the spatial structure of the mosaic into
% account.  The implementation does not make sense for stimuli that are not spatially
% uniform.
%
% 8/27/13  dhb  Pulled this out on its own.  Modularize, modularize, modularize.

%% Consistency checks
if (~strcmp(blankConeResponses.type,testConeResponses.type))
    error('Different type for passed blank and test structures.  Nope.');
end
if (length(blankConeResponses.coneNumbersToUse) ~= length(testConeResponses.coneNumbersToUse))
    error('Different length coneNumbersToUse fields passed blank and test structures.  Nope.');
end
if (any(blankConeResponses.coneNumbersToUse ~= testConeResponses.coneNumbersToUse))
    error('Different coneNumbersToUse fields passed blank and test structures.  Nope.');
end
if (length(blankConeResponses.oneConeEachClassStartIndices) ~= length(testConeResponses.oneConeEachClassStartIndices))
    error('Different length oneConeEachClassStartIndices fields passed blank and test structures.  Nope.');
end
if (any(blankConeResponses.oneConeEachClassStartIndices ~= testConeResponses.oneConeEachClassStartIndices))
    error('Different oneConeEachClassStartIndices fields passed blank and test structures.  Nope.');
end


%% Initialize output structures
blankSecondSiteResponses.type = blankConeResponses.type;
blankSecondSiteResponses.surroundType = theParams.surroundType;
blankSecondSiteResponses.theVectors = blankConeResponses.theVectors;
blankSecondSiteResponses.coneNumbersToUse = blankConeResponses.coneNumbersToUse;
blankSecondSiteResponses.oneConeEachClassStartIndices = blankConeResponses.oneConeEachClassStartIndices;

testSecondSiteResponses.type = testConeResponses.type;
testSecondSiteResponses.surroundType = theParams.surroundType;
testSecondSiteResponses.theVectors = testConeResponses.theVectors;
testSecondSiteResponses.coneNumbersToUse = testConeResponses.coneNumbersToUse;
testSecondSiteResponses.oneConeEachClassStartIndices = testConeResponses.oneConeEachClassStartIndices;

%% Do the second site calculations
switch(blankSecondSiteResponses.type)
    case 'rgb_uniform'
        switch (blankSecondSiteResponses.surroundType)
            case 'none'
                % No surround. Do nothing here.
                
            case 'random_wiring'
                % Subract a random mix of L and M cone responses
                % from each L and M cone response.
                %
                % This code doesn't exclude the center cone from the
                % possible cones from which the surround is computed.
                % I think this is OK for our special case of uniform fields,
                % even though it isn't exactly correct.
                %
                % The weight is applied to the mean response of the surround cones,
                % so that its magnitude remains commensurate across variation in the
                % number that are averaged.
                %
                % I learned the hard way that some care is needed to make sure this is fast.
                if (theParams.surroundSize > 0 && theParams.surroundWeight > 0)
                    centerIndices = 1:blankSecondSiteResponses.oneConeEachClassStartIndices(3)-1;
                    surroundIndices = 1:blankSecondSiteResponses.oneConeEachClassStartIndices(3)-1;
                    
                    % Each time through the surround size loop, radomize the order of the LM cones and
                    % subtract one weighted surround cone from the center.
                    %
                    % Blank.  Compute identity of each surround cone in this loop, and store.  Use
                    % same identities in the test loop below.
                    theOutVectors = blankSecondSiteResponses.theVectors;
                    for s = 1:theParams.surroundSize
                        % These two lines do the shuffle, and keep the dimensions matched for the subtract
                        % even when there are different numbers of center and surround cones.
                        otherConesIndicesIndices = Ranint(length(centerIndices),length(surroundIndices));
                        otherConesIndices{s} = surroundIndices(otherConesIndicesIndices);
                        theOutVectors(centerIndices,:) = theOutVectors(centerIndices,:) - ...
                            (theParams.surroundWeight/theParams.surroundSize)*blankSecondSiteResponses.theVectors(otherConesIndices{s},:);          
                    end
                    blankSecondSiteResponses.theVectors = theOutVectors;
                    
                    % Test
                    theOutVectors = testSecondSiteResponses.theVectors;
                    for s = 1:theParams.surroundSize
                        theOutVectors(centerIndices,:) = theOutVectors(centerIndices,:) - ...
                            (theParams.surroundWeight/theParams.surroundSize)*testSecondSiteResponses.theVectors(otherConesIndices{s},:);          
                    end
                    testSecondSiteResponses.theVectors = theOutVectors;
                end
                
            case 'cone_specific'
                % This subtracts the response of M cones from each L cone, and vice-versa.  Logic is the
                % same as above, but done twice, once for L centers and once for M centers.
                if (theParams.surroundWeight > 0)
                    
                    % L cone centers, blank
                    centerIndices = 1:blankSecondSiteResponses.oneConeEachClassStartIndices(2)-1;
                    surroundIndices = blankSecondSiteResponses.oneConeEachClassStartIndices(2):blankSecondSiteResponses.oneConeEachClassStartIndices(3)-1;
                    theOutVectors = blankSecondSiteResponses.theVectors;
                    for s = 1:theParams.surroundSize
                        otherConesIndicesIndices = Ranint(length(centerIndices),length(surroundIndices));
                        otherConesIndices{s} = surroundIndices(otherConesIndicesIndices);
                        theOutVectors(centerIndices,:) = theOutVectors(centerIndices,:) - ...
                            (theParams.surroundWeight/theParams.surroundSize)*blankSecondSiteResponses.theVectors(otherConesIndices{s},:);          
                    end
                    blankSecondSiteResponses.theVectors = theOutVectors;
                    
                    % L cone centers, test
                    theOutVectors = testSecondSiteResponses.theVectors;
                    for s = 1:theParams.surroundSize
                        theOutVectors(centerIndices,:) = theOutVectors(centerIndices,:) - ...
                            (theParams.surroundWeight/theParams.surroundSize)*testSecondSiteResponses.theVectors(otherConesIndices{s},:);          
                    end
                    testSecondSiteResponses.theVectors = theOutVectors;

                    % M cone centers, blank
                    centerIndices = blankSecondSiteResponses.oneConeEachClassStartIndices(2):blankSecondSiteResponses.oneConeEachClassStartIndices(3)-1;
                    surroundIndices = 1:blankSecondSiteResponses.oneConeEachClassStartIndices(2)-1;
                    theOutVectors = blankSecondSiteResponses.theVectors;
                    for s = 1:theParams.surroundSize
                        otherConesIndicesIndices = Ranint(length(centerIndices),length(surroundIndices));
                        otherConesIndices{s} = surroundIndices(otherConesIndicesIndices);
                        theOutVectors(centerIndices,:) = theOutVectors(centerIndices,:) - ...
                            (theParams.surroundWeight/theParams.surroundSize)*blankSecondSiteResponses.theVectors(otherConesIndices{s},:);
                    end
                    blankSecondSiteResponses.theVectors = theOutVectors;
                    
                    % M cone centers, test
                    theOutVectors = testSecondSiteResponses.theVectors;
                    for s = 1:theParams.surroundSize
                        theOutVectors(centerIndices,:) = theOutVectors(centerIndices,:) - ...
                            (theParams.surroundWeight/theParams.surroundSize)*testSecondSiteResponses.theVectors(otherConesIndices{s},:);
                    end
                    testSecondSiteResponses.theVectors = theOutVectors;
                end
                
            otherwise
                error('Unknown surround type specified for stimulus type ''rgb_uniform''');
        end
        
        %% Add noise at a second site
        %
        % The noise is governed by the specified Fano factor,
        % although we compute the variance from the magnitude
        % of the sample value, rather from a mean over time.  
        % That is mostly a matter of convenience, and perhaps
        % not unreasonable biophysically in any case (I'm not sure).
        if (theParams.secondSiteFanoFactor > 0)
            theVarMat = theParams.secondSiteFanoFactor*abs(blankSecondSiteResponses.theVectors);
            theSdMat = sqrt(theVarMat);
            blankSecondSiteResponses.theVectors = blankSecondSiteResponses.theVectors + theSdMat.*randn(size(blankSecondSiteResponses.theVectors));
            
            theVarMat = theParams.secondSiteFanoFactor*abs(testSecondSiteResponses.theVectors);
            theSdMat = sqrt(theVarMat);
            testSecondSiteResponses.theVectors = testSecondSiteResponses.theVectors + theSdMat.*randn(size(testSecondSiteResponses.theVectors));
        end
        
    otherwise
        error('Unknown cone response type passed');
end

