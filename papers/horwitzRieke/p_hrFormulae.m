%% Formulae from HR manuscript

%% Macaque
%  All cones per mm2, after Packer 1989 and Goodchild 1996
%  x in mm, we think
% 3 mm is
x     = 0:.05:3;
cpmm2 = 1e3*(150.9*exp(-1.2*x) + 35.9*exp(-0.15*x) + 9.9*exp(-0.03*x));

vcNewGraphWin;
plot(x,cpmm2); grid on

xlabel('position in mm');
ylabel('Cones per mmm^2');

%% In macaque mm / deg is 0.2
deg  = x/.2;
vcNewGraphWin;
plot(deg,cpmm2); grid on

xlabel('position in deg');
ylabel('Cones per mmm^2');

%% S-cones/mm2

scpmm2 = 2469.09*exp(-0.2*x) + 1822.54*exp(-0.05*x);
vcNewGraphWin;
plot(deg,scpmm2,'b-');
xlabel('position in deg');
ylabel('Cones per mmm^2');

semilogy(deg,cpmm2,'r-',deg,scpmm2,'b--')
xlabel('position in deg');
ylabel('Cones per mmm^2');

% Separate out the LM by subtracting
lmpmm2 = cpmm2 - scpmm2;
vcNewGraphWin;
semilogy(deg,lmpmm2,'r-',deg,scpmm2,'b--')
xlabel('position in deg');
ylabel('Cones per mmm^2');

% Cone separation at some position
pos = 30
nConesPerMM = sqrt(lmpmm2(pos))
coneSpacingUM = (1/nConesPerMM)*1e3


%% Background isomerizations per sec per cone
%  L - 7131, M - 6017, S 1973 (page 10 of their manuscript)
