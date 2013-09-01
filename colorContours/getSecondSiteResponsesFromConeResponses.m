function secondSiteResponses = getSecondSiteResponsesFromConeResponses(coneResponses,theParams,staticParams,runtimeParams,staticComptutedValues)
% theSecondSiteResponses = getSecondSiteResponsesFromConeResponses(coneResponses,theParams,staticParams,runtimeParams,staticComptutedValues)
%
% Get a set of second site responses from the cone responses.  How the second
% site responses are computed depends both on the type of the coneResponses
% structure, and on variables specified in the fields of theParams.
%
% 8/27/13  dhb  Pulled this out on its own.  Modularize, modularize, modularize.

%% Initialize output structure
secondSiteResponses.type = coneResponses.type;
secondSiteResponses.surroundType = theParams.surroundType;
secondSiteResponses.theVectors = coneResponses.theVectors;
secondSiteResponses.coneNumbersToUse = coneResponses.coneNumbersToUse;
secondSiteResponses.oneConeEachClassStartIndices = coneResponses.oneConeEachClassStartIndices;

%% Do the second site calculations
switch(secondSiteResponses.type)
    case 'rgb_uniform'
        switch (secondSiteResponses.surroundType)
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
                if (theParams.surroundSize > 0 && theParams.surroundWeight > 0)
                    for k = 1:staticParams.nDrawsPerTestStimulus
                        for j = 1:secondSiteResponses.oneConeEachClassStartIndices(3)-1  
                            centerConeVal = secondSiteResponses.theVectors(j,k);
                            otherConesIndex = Shuffle(1:secondSiteResponses.oneConeEachClassStartIndices(3)-1);
                            surroundConeVals = coneResponses.theVectors(otherConesIndex(1:theParams.surroundSize),k);
                            secondSiteResponses.theVectors(j,k) = centerConeVal - theParams.surroundWeight*sum(surroundConeVals)/theParams.surroundSize;
                        end
                    end
                end
                
            case 'cone_specific'
                % This subtracts the response M cones from each L cone, and vice-versa
                if (theParams.surroundWeight > 0)
                    for k = 1:staticParams.nDrawsPerTestStimulus
                        % Do each L cone
                        for j = 1:secondSiteResponses.oneConeEachClassStartIndices(2)-1
                            centerConeVal = secondSiteResponses.theVectors(j,k);
                            otherConesIndex = Shuffle(secondSiteResponses.oneConeEachClassStartIndices(2):secondSiteResponses.oneConeEachClassStartIndices(3)-1);
                            surroundConeVals = coneResponses.theVectors(otherConesIndex(1:theParams.surroundSize),k);
                            secondSiteResponses.theVectors(j,k) = centerConeVal - theParams.surroundWeight*sum(surroundConeVals)/theParams.surroundSize;
                        end
                        
                        % Do each M cone
                        for j = secondSiteResponses.oneConeEachClassStartIndices(2):secondSiteResponses.oneConeEachClassStartIndices(3)-1
                            centerConeVal = secondSiteResponses.theVectors(j,k);
                            otherConesIndex = Shuffle(1:secondSiteResponses.oneConeEachClassStartIndices(2)-1);
                            surroundConeVals = coneResponses.theVectors(otherConesIndex(1:theParams.surroundSize),k);
                            secondSiteResponses.theVectors(j,k) = centerConeVal - theParams.surroundWeight*sum(surroundConeVals)/theParams.surroundSize;
                        end
                    end
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
            theVarMat = theParams.secondSiteFanoFactor*abs(secondSiteResponses.theVectors);
            theSdMat = sqrt(theVarMat);
            secondSiteResponses.theVectors = secondSiteResponses.theVectors + theSdMat.*randn(size(secondSiteResponses.theVectors));
        end
        
    otherwise
        error('Unknown cone response type passed');
end

