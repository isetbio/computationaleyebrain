function LMScontMatrix = cm_ledamp2percont(Amps, LMSfunctions, meanLEDSPDs, meanLEDamps)
%
% LMScontMatrix = cm_ledamp2percont(Amps, LMSfunctions, meanLEDSPDs);
%
% cm_ledamp2percont is bacically working for 6LED experiments when you
% want to have LMS cone (with or without ipRGC) contrast responses.
% 
%<Input>
% Amps ... amplitudes of LEDs, which should be L x n matrix
%
% LMSfunctions ... LMS response functions, which should be w x s matrix.
%
% meanLEDSPDs  ... mean of LED spectrum power distribution, which should be w
%                  x L matrix.
%
% ,where L is number of color channel of display, n is number of amplitude
% set of display, w and s are wavelength and number of sensors,
% respectively.
%
%<Output>
% LMScontMatrix ... contrast of LMS responses, which should be s x n matrix
%
% % Exapmle
%
%

if ~exist('meanLEDSPDs','var'), meanLEDSPDs = cm_getledSPD(cm_getDefaultWls); end
if ~exist('meanLEDamps','var'), meanLEDamps = [.5 .5 .5 .5 .5 .5]; end

meanLEDamps = meanLEDamps * 2;

[numLEDs numPairs] = size(Amps);

meanSPD       = meanLEDSPDs   * (ones(numLEDs,1) .* meanLEDamps(:));   % meanSPD is w x 1 vector
meanResponse  = LMSfunctions' * meanSPD;           % meanResponse is s x 1 vector
meanResponses = meanResponse  * ones(1, numPairs); % meanResponses is s x n matrix

SPDs          = meanLEDSPDs   * Amps;              % SPD is w x n matrix
% Column of SPDs corresponds with Column of Amps/  

LMSresponse   = LMSfunctions' * SPDs;              % LMSresponse is s x n matrix
% Column of LMSresponse also corresponds with Column of Amps 

LMScontMatrix = LMSresponse  ./ meanResponses;

end
