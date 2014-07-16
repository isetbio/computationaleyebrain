function classificationData = oneIntervalToTwoIntervalClassificationData(classificationData)
% classificationData = oneIntervalToTwoIntervalClassificationData(classificationData)
%
% Use the data in the classification structure, set up for a one-interval experiment and
% convert it into a form for TAFC.
%
% 8/29/13  dhb  Pulled this out into its own routine.

% Build up a TAFC training and test set from the one interval data
oneIntervalDataDimension = size(classificationData.trainingData,2);
nTrainingData = 2*size(classificationData.trainingData,1);
trainingBlankIndices = find(classificationData.trainingLabels == classificationData.blankLabel);
trainingTestIndices = find(classificationData.trainingLabels == classificationData.testLabel);
validateBlankIndices = find(classificationData.validateLabels == classificationData.blankLabel);
validateTestIndices = find(classificationData.validateLabels == classificationData.testLabel);
tafcTrainingData = zeros(nTrainingData,2*oneIntervalDataDimension);
tafcTrainingLabels = zeros(nTrainingData,1);
tafcValidateData = zeros(nTrainingData,2*oneIntervalDataDimension);
tafcValidateLabels = zeros(nTrainingData,1);
for tt = 1:nTrainingData
    % Training set
    %
    % Flip a coin to decide whether test is in first or second interval
    if (CoinFlip(1,0.5))
        temp1 = Shuffle(trainingTestIndices);
        temp2 = Shuffle(trainingBlankIndices);
        tafcTrainingLabels(tt) = classificationData.blankLabel;
    else
        temp1 = Shuffle(trainingBlankIndices);
        temp2 = Shuffle(trainingTestIndices);
        tafcTrainingLabels(tt) = classificationData.testLabel;
    end
    tafcTrainingData(tt,1:oneIntervalDataDimension) = classificationData.trainingData(temp1(1),:);
    tafcTrainingData(tt,oneIntervalDataDimension+1:2*oneIntervalDataDimension) = classificationData.trainingData(temp2(1),:);
    
    % Validation set
    %
    % Same logic
    if (CoinFlip(1,0.5))
        temp1 = Shuffle(validateTestIndices);
        temp2 = Shuffle(validateBlankIndices);
        tafcValidateLabels(tt) = classificationData.blankLabel;
    else
        temp1 = Shuffle(validateBlankIndices);
        temp2 = Shuffle(validateTestIndices);
        tafcValidateLabels(tt) = classificationData.testLabel;
    end
    tafcValidateData(tt,1:oneIntervalDataDimension) = classificationData.validateData(temp1(1),:);
    tafcValidateData(tt,oneIntervalDataDimension+1:2*oneIntervalDataDimension) = classificationData.validateData(temp2(1),:);
end

%% Rewrite the structure fields for return
classificationData.trainingData = tafcTrainingData;
classificationData.trainingLabels = tafcTrainingLabels;
classificationData.validateData = tafcValidateData;
classificationData.validateLabels = tafcValidateLabels;

