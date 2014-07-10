function varLMSIfun = cm_variableLMSI_PODandLambda(AbsSpectra, POD, Lambdashift, eyetransmittance)
%
% varLMSIfun = cm_variableLMSI_PODandLambda(AbsSpectra, POD, Lambdashift, eyetransmittance)
% 
% <Input>
%
% AbsSpectra        ... w x n, absorbance spectrum power distribution
% POD               ... 1 x n, photopigment optical density
% Lambdashift       ... 1 x n, Lambda max shift (nm)
% eyetransmittance  ... 1 x n, whole eye transmittance (macular and lens)
%
% <Output>
%
% varLMSIfun        ... w x n, LMSI response function based on axial photopigment density and shift of lambda max 
%
%
% C) VistaLab 2012 HH 
%
%%
if ~exist('wls','var'),              wls = cm_getDefaultWls;                     end
if ~exist('eyetransmittance','var'), eyetransmittance = LensTransmittance(wls'); end
    
if length(POD) == 1;            POD = POD * ones(1,4);                  end

if ~exist('Lambdashift','var') || isempty(Lambdashift)    
       skipshiftflag = true;
else
       skipshiftflag = false;
       if length(Lambdashift) == 1;
           Lambdashift = Lambdashift * ones(1,4);
           skipshiftflag = false;
       elseif Lambdashift == 0, skipshiftflag = true;          
       end   
end

numSensor = size(AbsSpectra,2);

for ii = 1:numSensor
    % Absorbance to Absorbtance
    varAbtLMSI = AbsorbanceToAbsorbtance(AbsSpectra(:,ii)', wls, POD(ii));
    
    if skipshiftflag == false
        % shift lambda max of Absorbtance
        varAbtLMSI = interp1(wls + Lambdashift(ii) ,varAbtLMSI, wls, 'cubic', 0); % use interp
    end
    
    % Absorbtance to Response function
    tmp = varAbtLMSI' .* eyetransmittance(:);
        
    % normalize (max relative response should be one.)
    varLMSIfun(:,ii) = tmp ./ max(tmp(:));
    
end


end

