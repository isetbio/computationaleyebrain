
%% Brettel appendix calculation of the cones

% The way Brettel and people calculate the cones is kind of nuts because
% they include scale factors that depend on the display monitor.  Still, we
% follow John Mollon like the religious figure that he is.

% Find k depending on the white point of the monitor
w = 400:700;
d = displayCreate('CRT-Dell',400:700);
whiteSPD = displayGet(d,'white spd');

% Read the Stockman data
stockman = ieReadSpectra('stockman',w);
stockman = stockman*diag([0.68273, 0.35235, 1]);

% Here is how they scale.  This is intended to make L+M kind of like
% Vlambda.  See Equation 1
n1 = sum(stockman*[1 1 0]');
n2 = sum(stockman(:,3));
stockman = stockman*(diag(1./[n1 n1 n2]));
vcNewGraphWin; plot(w,stockman)

% Here is the magic scaling by the display
whiteLMS = whiteSPD'*stockman
k = whiteLMS(1) + whiteLMS(2);
stockman = stockman / k;

whiteSPD'*stockman

%% Here are the anchor points
pts = [475 485 476 660];

for ii = pts
    idx = find(w == ii)
    stockman(idx,:)
end

