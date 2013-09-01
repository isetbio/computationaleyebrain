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
                %
                % I learned the hard way that some care is needed to make sure this is fast.
                if (theParams.surroundSize > 0 && theParams.surroundWeight > 0)
                    centerIndices = 1:secondSiteResponses.oneConeEachClassStartIndices(3)-1;
                    surroundIndices = 1:secondSiteResponses.oneConeEachClassStartIndices(3)-1;
                    
                    % Each time through the surround size loop, radomize the order of the LM cones and
                    % subtract one weighted surround cone from the center.
                    for s = 1:theParams.surroundSize
                        % These two lines do the shuffle, and keep the dimensions matched for the subtract
                        % even when there are differnt numbers of center and surround cones.
                        otherConesIndicesIndices = Ranint(length(centerIndices),length(surroundIndices));
                        otherConesIndices = surroundIndices(otherConesIndicesIndices);
                        secondSiteResponses.theVectors(centerIndices,:) = secondSiteResponses.theVectors(centerIndices,:) - ...
                            (theParams.surroundWeight/theParams.surroundSize)*secondSiteResponses.theVectors(otherConesIndices,:);          
                    end
                end
                
            case 'cone_specific'
                % This subtracts the response of M cones from each L cone, and vice-versa.  Logic is the
                % same as above, but done twice, once for L centers and once for M centers.
                if (theParams.surroundWeight > 0)
                    
                    % L cone centers
                    centerIndices = 1:secondSiteResponses.oneConeEachClassStartIndices(2)-1;
                    surroundIndices = secondSiteResponses.oneConeEachClassStartIndices(2):secondSiteResponses.oneConeEachClassStartIndices(3)-1;
                    for s = 1:theParams.surroundSize
                        otherConesIndicesIndices = Ranint(length(centerIndices),length(surroundIndices));
                        otherConesIndices = surroundIndices(otherConesIndicesIndices);
                        secondSiteResponses.theVectors(centerIndices,:) = secondSiteResponses.theVectors(centerIndices,:) - ...
                            (theParams.surroundWeight/theParams.surroundSize)*secondSiteResponses.theVectors(otherConesIndices,:);          
                    end
                    
                    % M cone centers
                    centerIndices = secondSiteResponses.oneConeEachClassStartIndices(2):secondSiteResponses.oneConeEachClassStartIndices(3)-1;
                    surroundIndices = 1:secondSiteResponses.oneConeEachClassStartIndices(2)-1;
                    for s = 1:theParams.surroundSize
                        otherConesIndicesIndices = Ranint(length(centerIndices),length(surroundIndices));
                        otherConesIndices = surroundIndices(otherConesIndicesIndices);
                        secondSiteResponses.theVectors(centerIndices,:) = secondSiteResponses.theVectors(centerIndices,:) - ...
                            (theParams.surroundWeight/theParams.surroundSize)*secondSiteResponses.theVectors(otherConesIndices,:);          
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

