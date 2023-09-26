addpath('/Users/hopelutwak/Documents/MATLAB/VisTools/')
    view_dist                       = .35;                                  % m .57; psychophysics room: .35 .50
    screensize                      = [.405 .298];                           % m
    pixels                          = [1920 1200];
    frame_rate                      = 60;   
    stimulus_duration               = 1;
    
    view_window                     = round([rad2deg(atan(screensize(1)/2/view_dist)) rad2deg(atan(screensize(2)/2/view_dist))]);  % X,Y centered around fixation [60 46] [54 40];                                                                      % [36 27] for laptop, [55 32] for monitor
    scale_factor                    = sqrt((view_window(1)*60*view_window(2)*60)/(pixels(1)*pixels(2)))*1;          % Arcmin/pixel 1.78 2.25 1.92dell: 1680x1050, psychophysics room: 1600x1200                                % Screen frame rate (hz) psychophysics room: 85
    scale_factorX = view_window(1)*60/pixels(1);
    scale_factorY = view_window(2)*60/pixels(2);
   
    
    text_diam                         = 60; 
    text_pixelW                       = ceil(text_diam/scale_factor);

    H_ecc_fix                       = 0;                                    % Horizontal fixation ecc (degs, neg is left)
    V_ecc_fix                       = 0;                                    % Vertical fixation ecc (degs, neg is up)

    H_ecc_fix                       = H_ecc_fix*60/scale_factorX;
    V_ecc_fix                       = V_ecc_fix*60/scale_factorY;
    mv_length                       = ceil(stimulus_duration*frame_rate);

    bps = text_diam+1;
    numSame                         = 1;                                   % number of apertures that have same texture
    cut                           = 3;                                    %cut off first 3rd frequencies
    
%% OPTION 1: drifting texture 
stimulus_duration = 1;
frame_rate = 60;

velocity_field = [[1,1];[1,1]]; %deg/s

dots_deg = [5, -5; 0, 0];

%%
%pixels per s
% velocity_field_pix = [velocity_field(1,:)*60*1/scale_factorX; velocity_field(2,:)*60*1/scale_factorX]; %convert to pixel/sec
velocity_field_pix = velocity_field*60*1/scale_factor;
%convert to pixel/frame
velocity_field_ppf = velocity_field_pix/frame_rate;

% sz = ceil(stimulus_duration*frame_rate) * max(speeds) + text_pixelW;
% sz = max(speeds)*stimulus_duration + text_pixelW;
sz = max(vecnorm(velocity_field_pix))*2*stimulus_duration + text_pixelW;
numDrifts = ceil(length(velocity_field)/numSame);

%%
Screen('Preference', 'SkipSyncTests', 1);
AssertOpenGL;
screens=Screen('Screens');
screenNumber=max(screens);

w=Screen('OpenWindow',screenNumber, 128);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
screen_rect = Screen('Rect',w);

%-----Screen Landmarks------------
sr_hor = round(screen_rect(3)/2); % Middle of the screen, horizontally, in pixels
sr_ver = round(screen_rect(4)/2); % Middle of the screen, vertically, in pixels
fix_hor = sr_hor+H_ecc_fix;     % Horizontal location of fixation cross, in pixels
fix_ver = sr_ver+V_ecc_fix;     % Vertical location of fixation cross, in pixels
movie_rect= [0,0,bps,bps];



tex = zeros(1,numDrifts);
for d = 1:numDrifts
    drift = oneoverfcut(1, ceil(sz),ceil(sz),cut);
    tex(d)=Screen('MakeTexture', w, 255*drift);
end
drift_center = [size(drift,1)/2, size(drift,1)/2];


% windowTex = Screen('MakeTexture',w,cat(3,128*ones(2*W),255*(1-ggaus(2*W,.2*W))));
windowTex = Screen('MakeTexture',w,cat(3,128*ones(text_pixelW),255*(1-coswin(text_pixelW,.3*text_pixelW,.5*text_pixelW))));

pos0 = CenterRectOnPoint(movie_rect, drift_center(1), drift_center(1));
pos = CenterRectOnPoint(movie_rect,dots_deg(1,:)'*60/scale_factorX+sr_hor,dots_deg(2,:)'*60/scale_factorY+sr_ver)'; %big Y numbers are lower on the screen in pixels so we have to flip

% source = cell(1,length(stimulus_duration*frame_rate));
% % calculate source rects for different speeds
% 
% for tt = 1:stimulus_duration*frame_rate
%     for s = 1:length(velocity_field)
%         source{tt}(:,s) = ...
%             [drift_center(1)-text_pixelW/2-(tt-1)*velocity_field_ppf(1,s),drift_center(2)-(tt-1)*velocity_field_ppf(2,s), ...
%             drift_center(1)+text_pixelW/2-(tt-1)*velocity_field_ppf(1,s),drift_center(2)+text_pixelW/2-(tt-1)*velocity_field_ppf(2,s)]';
%     end
% end
source = cell(1,length(stimulus_duration*frame_rate));
for tt = 1:stimulus_duration*frame_rate
    for s = 1:length(velocity_field)
        source{tt}(:,s) = [pos0(1)-(tt-1)*velocity_field_ppf(1,s); pos0(2)-(tt-1)*velocity_field_ppf(2,s); ...
            pos0(3)-(tt-1)*velocity_field_ppf(1,s); pos0(4)-(tt-1)*velocity_field_ppf(2,s)];
    end
end
% play movie
which_tex = randi(numDrifts, 1,length(velocity_field)); %assign one of drifts to each patch for this movie

for tt=1:stimulus_duration*frame_rate
    Screen('DrawTextures', w, tex(which_tex), source{tt} ,pos ,0);
    Screen('DrawTextures', w, windowTex, [], pos,0);
    Screen('Flip', w);
end

Screen('Close');
sca;


%% OPTION 2: precalculated movie

% clear;
% sca;
speed = 90; % pixels/s
stimulus_duration = 3;
frame_rate = 85;
text_pixelW = 120;

drift = 255*generateColoredNoiseMotion(text_pixelW,stimulus_duration,frame_rate,[speed,0],[0,0]);

Screen('Preference', 'SkipSyncTests', 1);
AssertOpenGL;
screens=Screen('Screens');
screenNumber=max(screens);

w=Screen('OpenWindow',screenNumber, 128);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % Turns on blending so you can use alpha values on images


for tt=1:stimulus_duration*frame_rate
   tex(tt) = Screen('MakeTexture', w, drift(:,:,tt)); 
end
windowTex = Screen('MakeTexture',w,cat(3,128*ones(2*text_pixelW),255*(1-ggaus(2*text_pixelW,.5*text_pixelW))));
% windowTex = Screen('MakeTexture',w,cat(3,128*ones(2*W),255*(1-coswin(2*W,.2*W,W))));
% windowTex = Screen('MakeTexture',w,cat(3,128*ones(W),255-circle));
angle = ones(size(angle))*45;
positions = [sr_hor+pos(1,:); sr_ver+pos(2,:); sr_hor+pos(3,:); sr_ver+pos(4,:)];
for tt=1:stimulus_duration*frame_rate    
    Screen('DrawTextures', w, tex(tt),[],pos, angle);
    Screen('DrawTextures', w, windowTex, [], pos, angle);
    Screen('Flip', w);
end

Screen('Close');
sca;
