function acc = svmClassifyAcc(dataMatrix, labels, nFolds, svmType, opts)
%% function getSVMAccuracy(dataMatrix, labels, [nFolds], [svmType], [opts])
%    compute the svm prediction accuracy
%
%  Required packages:
%    libsvm: free copy can be downloaded from
%            http://www.csie.ntu.edu.tw/~cjlin/libsvm/
%    liblinear: free copy can be downloaded from
%            http://www.csie.ntu.edu.tw/~cjlin/liblinear/
%
%  Inputs:
%    dataMatrix - M-by-N matrix, containing data for M instances and N
%                 features
%    labels     - M-by-1 vector, label for each instance
%    nFolds     - number of folds for cross validation
%    svmType    - string, can be 'svm' or 'linear', indicating whether to 
%                 use libsvm or liblinear for classification, by default,
%                 'svm' is used
%    opts       - string, libsvm / liblinear options
%
%  Outputs      - 2-by-1 vector, containing average accuray and standard
%                 deviation
%
%  Example:
%    acc = svmClassifyAcc(dataMatrix, labels, 10, 'linear')
%
%  See also:
%    svmtrain, svmpredict, train, predict
%
%  (HJ) VISTASOFT Team 2013

%% Check Inputs
if nargin < 1, error('Data matrix required'); end
if nargin < 2, error('Labels for data required'); end
if nargin < 3, nFolds  = 5; end
if nargin < 4, svmType = 'svm'; end
if nargin < 5, opts = '-s 2 -q'; end
if isempty(opts), opts = '-s 2 -q'; end

%% Divide data into nFolds
%  Normalize Data
%  Normalization is important when data values are small, i.e. volts image
%  If we values are large, i.e. photon absorptions, it's fine to skip this
%  step
dataMatrix = (dataMatrix-repmat(min(dataMatrix),[size(dataMatrix, 1) 1])) ...
    ./ repmat(max(dataMatrix)-min(dataMatrix),[size(dataMatrix,1 ) 1]);

% The data should not contain any NaN inside
% If there's any NaN, set it to zero
dataMatrix(isnan(dataMatrix)) = 0;

%  Random permute the data
[M, ~] = size(dataMatrix);
ind   = randperm(M);

% Train and Test
accHistory  = zeros(nFolds,1);
instPerFold = round(M/nFolds);
for i = 1 : nFolds
    if i < nFolds
        trainIndx = [ind(1:(i-1)*instPerFold) ...
                     ind(i*instPerFold+1:end)];
        testIndx  = ind((i-1)*instPerFold+1:i*instPerFold);
    else
        trainIndx = ind(1:(i-1)*instPerFold);
        testIndx  = ind((i-1)*instPerFold+1:end);
    end
    trainData = sparse(dataMatrix(trainIndx,:));
    testData  = sparse(dataMatrix(testIndx,:));
    % Train
    switch svmType
        case 'linear' 
        % Liblinear routine
          svmStruct = train(labels(trainIndx),trainData,opts);
          [~,curAcc,~] = predict(labels(testIndx),testData,svmStruct,'-q');
        case 'svm'
        % LibSVM routine
        % Parameters explaination:
        %   http://www.csie.ntu.edu.tw/~cjlin/libsvm/
          svmStruct = svmtrain(labels(trainIndx),trainData,opts);
          [~,curAcc,~] = svmpredict(labels(testIndx), testData, ...
                                    svmStruct,'-q');
        otherwise
          error('Unknown svm type');
    end
    accHistory(i) = curAcc(1) / 100; % Convert to between 0~1
end

% Report average and std
acc(1) = mean(accHistory);
acc(2) = std(accHistory);


end