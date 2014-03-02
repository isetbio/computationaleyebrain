
%% Brettel appendix calculation of the cones

w = 400:700;
stockman = ieReadSpectra('stockman',w);
stockman = stockman*diag([0.68273, 0.35235, 1]);

n1 = sum(stockman*[1 1 0]');
n2 = sum(stockman(:,3));

stockman = stockman*(diag(1./[n1 n1 n2]));
vcNewGraphWin; plot(w,stockman)

% Find k depending on the white point of the monitor


ii = find(w == 485)
stockman(ii,:)
