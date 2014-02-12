function classificationData = getTrainingAndValidationData(blankResponses,testResponses)
% classificationData = getTrainingAndValidationData(blankResponses,testResponses)
%
% Get training and validation data sets out of the test and blank
% responses, and put into a structure for use in the classification.
%
% [**] It would be good to document the structure.
%
% 8/27/13  dhb  Wrote it.

%% Get and check consistency of stimulus type
if (~strcmp(blankResponses.type,testResponses.type))
    error('Blank and test not generated from consistent stimulus types');
end

%% Process according to type
switch (blankResponses.type)
    case 'rgb_uniform'
        %% Get number of cone classes
        nConeClasses = length(blankResponses.coneNumbersToUse);
        if (length(testResponses.coneNumbersToUse) ~= nConeClasses)
            error('Inconsistent length of coneNumbersToUse field across blank and test');
        end
        
        %% Get starting indices for each cone class, and consistency checks
        oneConeEachClassStartIndices = blankResponses.oneConeEachClassStartIndices;
        if (length(oneConeEachClassStartIndices) ~= nConeClasses)
            error('Number of cone classes not consistent across places where it matters.');
        end
        if (any(oneConeEachClassStartIndices ~= blankResponses.oneConeEachClassStartIndices))
            error('Cone class start indices not consistent across blank and test');
        end
        
        %% Get number of draws, number of cones, and consistency checks
        coneNumbersToUse = blankResponses.coneNumbersToUse;
        if (any(testResponses.coneNumbersToUse ~= coneNumbersToUse))
            error('Number of cones from each class to use not consistent across blank and test');
        end
        numberConesPerDraw = sum(coneNumbersToUse);
        
        numberOfCones = size(blankResponses.theVectors,1);
        if (size(testResponses.theVectors,1) ~= numberOfCones)
            error('Number of cones simulated differs between blank and test');
        end
        numberDraws = size(blankResponses.theVectors,2);
        if (size(testResponses.theVectors,2) ~= numberDraws)
            error('Blank and test data do not have same number of draws.');
        end
        
        %% Fill in the data
        %
        % Set up start indices
        startInIndex = oneConeEachClassStartIndices(1);
        if (startInIndex ~= 1)
            error('First class start index is not one.  This might be right.  But it is surprising and might indicated a bug.');
        end
        startOutIndex = 1;
        
        % Allocate space
        blankData = zeros(numberConesPerDraw,numberDraws);
        testData = zeros(numberConesPerDraw,numberDraws);
        
        % Select the right amount of data
        for c = 1:nConeClasses
            % End index for input
            if (c < nConeClasses)
                endInIndex = oneConeEachClassStartIndices(c+1) - 1;
            else
                endInIndex = numberOfCones;
            end
            
            % Randomly select the right number of cones for this class
            inIndicesToUse = Shuffle(startInIndex:endInIndex);
            inIndicesToUse = inIndicesToUse(1:coneNumbersToUse(c));
            
            % End index for output
            endOutIndex = startOutIndex+coneNumbersToUse(c)-1;
            blankData(startOutIndex:endOutIndex,:) = blankResponses.theVectors(inIndicesToUse,:);
            testData(startOutIndex:endOutIndex,:) = testResponses.theVectors(inIndicesToUse,:);
            
            % Bump start indices
            startInIndex = endInIndex+1;
            startOutIndex = endOutIndex+1;
        end
        
        
        %% Pull apart into training and validation datasets and
        % put into a structure for classification use.
        % 
        % As part of this, we transpose the data set, so that each row
        % represents a draw, with the data across the columns.  This is
        % the form that svm wants, so we might as well do it here.
        classificationData.blankLabel = -1;
        classificationData.testLabel = 1;
        fullData = [blankData' ; testData'];
        fullDataN = size(fullData,1);
        trainingDataN = round(fullDataN/2);
        validateDataN = fullDataN-trainingDataN;
        fullLabels = [classificationData.blankLabel*ones(trainingDataN,1) ; classificationData.testLabel*ones(validateDataN,1)];
        indices = Shuffle(1:size(fullData,1));
        classificationData.trainingData = fullData(indices(1:trainingDataN),:);
        classificationData.trainingLabels = fullLabels(indices(1:trainingDataN));
        classificationData.validateData = fullData(indices(trainingDataN+1:end),:);
        classificationData.validateLabels = fullLabels(indices(trainingDataN+1:end));
        
    otherwise
        error('Unknown stimulus type specified');
end