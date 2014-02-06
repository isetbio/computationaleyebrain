%% s_emFitParameters
%    This script is used to analyze the eye movement data from Jonathan
%    Winawerb.
%
%    In this script, we will fit a Bronian motion model with bounce back
%
%  (HJ) Jan, 2014

%% Clean up
clear all; close all; clc;

%% Load Data
%  Load original data
emData = importdata('jon1_samples.txt');

%  Get X, Y position
xPos = emData(:,2); yPos = emData(:,3);

%  Get rid of NAN
indx = ~isnan(xPos) & ~isnan(yPos);
xPos = xPos(indx); yPos = yPos(indx);

%  Get rid of extreme data by using 1% to 99% quantile
xMax = quantile(xPos, 0.99); xMin = quantile(xPos, 0.01);
yMax = quantile(yPos, 0.99); yMin = quantile(yPos, 0.01);
indx = (xPos > xMin) & (xPos < xMax) & (yPos > yMin) & (yPos < yMax);

xPos = xPos(indx); yPos = yPos(indx);

%  Convert to degrees
xPos = (xPos - 640) * 0.038; yPos = (yPos - 512) * 0.038;

%  Plot
vcNewGraphWin;
plot(xPos); hold on; plot(yPos, 'r');

%% Fit Brownian motion
xDiff = xPos(2:end) - xPos(1:end-1);
yDiff = yPos(2:end) - yPos(1:end-1);

%  Get rid of non-continuity
%  Since we delete some extreme data and NaNs, there might be some
%  discontinuity in xPos and yPos. This will lead to extreme value in xDiff
%  and yDiff
xMax = quantile(xDiff, 0.99); xMin = quantile(xDiff, 0.01);
yMax = quantile(yDiff, 0.99); yMin = quantile(yDiff, 0.01);
indx = (xDiff > xMin) & (xDiff < xMax) & (yDiff > yMin) & (yDiff < yMax);

xDiff = xDiff(indx); yDiff = yDiff(indx);

%  Compute mean and covariance matrix
mu = [mean(xDiff) mean(yDiff)];
Sigma = cov([xDiff yDiff]);

disp('mean:'); disp(mu);
disp('covariance:'); disp(Sigma);

%% Evaluate how well the data fits Gaussian
skewness([xDiff yDiff])
kurtosis([xDiff yDiff])
