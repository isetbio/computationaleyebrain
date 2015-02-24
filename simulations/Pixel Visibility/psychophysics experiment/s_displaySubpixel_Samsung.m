%% s_displaySubpixel_Samsung
%
%
%  (HJ) ISETBIO TEAM, 2015

%% Create display structure
d = displayCreate('LCD-Gechic');
% vcAddObject(d); displayWindow;

%% Init display optional parameters
d.frameRate = 60; % refresh rate of the display

d.screenNumber = max(Screen('Screens')); % get screen number

% set display resolution
% screenRes = Screen('Resolution', d.screenNumber);
% d.numPixels = [screenRes.width screenRes.height];
d.numPixels = [1920 1080];
d.backColorRgb = 255 * [0.5 0.5 0.5]'; % set background color to gray
d.numBuffers = 2;

%% Initialize parameters for experiments
stimParams = initStimParams;
% This is problematic because we don't use bit++ here
% should fix it 
d = initFixParams(d, 0.25); % fixation fov to 0.25
stairParams = initStaircaseParams;

% init feedback sound
correctSnd = soundFreqSweep(200, 500, .1);
incorrectSnd = soundFreqSweep(1000, 200, .5);

%% Custumize the instruction Page
instructions{1} = 'Pixelation Visibility Test\n';
instructions{2} = 'This test is composed of several trials. In each trial, you will be presented two pathes, one uniform and one pixelated\n';
instructions{3} = 'Your task is to tell which patch is uniform. Press Z for the left patch and X for the right one\n';
instructions{4} = 'Press any key to continue';

%% Open screen
%  open screen with desired resolution and framerate
d = openScreen(d, 'hideCursor', false);

%% Do experiments
%  Well, I give up using staircase controll here, because the stairs cannot
%  be well-defined in this case. The current setting is to loop over all
%  possible cases
%
%  Init parameters
quitKey = 'q';
device  = getBestDevice(d);

% Generate keyList for checking responses after the trial
% If you need to use some special keys, including number keys, use cell
% For example, {'1!','2@'}
keyList     = zeros(1,256);
includeKeys = zeros(1,length(stairParams.responseSet));
for i=1:size(stairParams.responseSet,2)
    if(iscell(stairParams.responseSet))
        includeKeys(i) = KbName(stairParams.responseSet{i});
    else
        includeKeys(i) = KbName(stairParams.responseSet(i));
    end
end

includeKeys = [includeKeys KbName(quitKey)]; % add quitKey
includeKeys = unique(includeKeys);

keyList(includeKeys) = 1;

% Show instructions to user
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
pressKey2Begin(d, 0, false, cell2mat(instructions));

% test colors
c = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1; ...
    .34 0 0; 0 .11 0; 0 0 1;];
c = reshape(c, [size(c,1), 1, size(c,2)]);

% patterns
pixelSz = 2:3;
nTrials = 25; % repeat per case

% sequence of experiment
idx = randperm(nTrials * length(pixelSz) * size(c, 1));

% init result matrix
result = nan(size(c,1), length(pixelSz), nTrials);

% start loop over cases
for ii = 1 : length(idx)
    [ic, ip, in] = ind2sub([size(c,1), length(pixelSz), nTrials], idx(ii));
    
    % create image
    I = zeros(d.resolution(2), d.resolution(1), 3);
    pos = rand > 0.5;
    if pos, dst = [201 201]; alt_dst = [201 1001]; 
    else dst = [201, 1001]; alt_dst = [201 201]; end
    kernel = zeros(2*pixelSz(ip), 2*pixelSz(ip), 3);
    kernel(1:pixelSz(ip), 1:pixelSz(ip),:) = ...
        repmat(c(ic,:,:), [pixelSz(ip) pixelSz(ip)]);
    I(dst(1):dst(1)+599,dst(2):dst(2)+599, :) = ...
        repmat(kernel, [600/pixelSz(ip)/2 600/pixelSz(ip)/2]);
    I(alt_dst(1):alt_dst(1)+599, alt_dst(2):alt_dst(2)+599, :) = ...
        repmat(c(ic,:,:), [600 600])/4;
    
    % create texture and show
    texturePtr = Screen('MakeTexture', d.windowPtr, I, [], [], 2);
    Screen('DrawTexture', d.windowPtr, texturePtr);
    Screen('Flip', d.windowPtr);
    
    % wait for response
    [~,KeyCode] = KbWait(device); WaitSecs(0.1);
    if KbName(KeyCode) == 'z', res = 1; else res = 0; end
    
    % record
    result(ic, ip, in) = (res == pos);
    
    % feedback
    if res == pos
        sound(correctSnd);
    else
        sound(incorrectSnd);
    end
    
    % draw black and wait
    Screen('FillRect', d.windowPtr, [0 0 0]);
    Screen('Flip', d.windowPtr);
    WaitSecs(0.5);
end                                                                                 

%% Close window and clean up
closeScreen(d);

%% Analyze and plot

