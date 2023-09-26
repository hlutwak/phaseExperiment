%% Phase experiment, but staircasing with linear changes along the constraint
% updated HL 9/6/2020
clear all;
close all;
clc;
% addpath(genpath('/Applications/Psychtoolbox'))
% addpath('/Users/hopelutwak/Documents/MATLAB/VisTools/')
% addpath('/Users/hopelutwak/Documents/opticflow_old/objectDetection/RangeDatabase1080p')

%vvv vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    %-----Subject Settings------------------ 
    initials                        = 'FH';                                 % Subject initials used to save file
    session                         = 'C1';
    trainingMode                    = 1;                                    % 1 = gives feedback after each trial (training), 0 = no feedback (actual experiment)

    %-----Rig Settings----------------------
    mouseclick                      = 1;                                    % use mouseclick instead of keys
    view_dist                       = .35;                                  % m .57; psychophysics room: .35 .50
    screensize                      = [.405 .298];                           % m
    pixels                          = [1600 1200];                          %home monitor 1920x1080 psychophysics[1600 1200]
    frame_rate                      = 30;                                   % psychophysics room 85

    view_window                     = round([rad2deg(atan(screensize(1)/2/view_dist)) rad2deg(atan(screensize(2)/2/view_dist))]);  % X,Y centered around fixation [60 46] [54 40];                                                                      % [36 27] for laptop, [55 32] for monitor
    scale_factor                    = sqrt((view_window(1)*60*view_window(2)*60)/(pixels(1)*pixels(2)))*1;          % Arcmin/pixel 1.78 2.25 1.92dell: 1680x1050, psychophysics room: 1600x1200                                % Screen frame rate (hz) psychophysics room: 85
    scale_factorX = view_window(1)*60/pixels(1);
    scale_factorY = view_window(2)*60/pixels(2);
    linearize                       = 0;                                    % Use calibrated LUT (do this when available)                    

    %-----Trial Settings------------------ 
    subjectNumber                   = '01';                                 % Number used to save eyelink file
    n_trials                        = 20;                                  % Number of trials
    translation                     = [[0, 0.15, 0.65]; [ 0, -0.15, .65]];             
%                                       [ 0, -0.15, .65]];                 % deg/s *up is down and down is up
    deviation_scale                 = 3;
    scramble                        = 1;                                  % scramble surround with optic flow info, 0 is same velociy in surround
    depth_structure                 = 0;                                    % 0 = random cloud, 1 = ground plane, 2 = natural scene 
    if depth_structure              ==1
        gaze_angle                  = 77;                                   %deg from staring straight at ground
        translation                 = getT_eyecoords(90-gaze_angle, translation);
        height                      = 1.65;                                 % height in meters
        bumpiness                   = 0.025;                                  % random noise in ground plane
        wall                        = 0;                                    % ground plane + wall starting in upper visual field
    elseif depth_structure          ==2
        load ('/Users/hopelutwak/Documents/opticflow_old/objectDetection/RangeDatabase1080p/lRange009')
    end
    
    theta                           = [0 45 90 135 180 225 270 315];        % range directions to test around base velocity
    fixateMovement                  = 0;                                    % 1 = red fixation dot moves, 0 = stationary red fixation dot
    eyeTracking                     = 0;                                    % 1 = eye tracking on, 0 = eye tracking off
    stimulus_duration               = .5;                                   % seconds
    test                            = 1;                                    % just shows a few stimv
    if test
        stimulus_duration = 1;
    end
    n_staircases                    = 2;                                    % Want to run multiple staircases
    linearize                       = 0;                                    % Use calibrated LUT (do this when available)
    data_path                       = './data/';
    
    % % Decides if experiment is simulated rotation or eye tracking and sets variables accordingly
    % linear_velocity                 = (rotation * .0174533) * 30;         % m/s, calculated from rotation (angular velocity)
    % if fixateMovement == 1
    %     rotation = 0;vv
    % elseif fixateMovement == 0
    %     linear_velocity = 0;
    % end
    %-----Experiment Settings, Don't Change---
    %-----Array Settings
    z0                              = 3;                                    % fixation distance
    middle_dist                     = 3; 
    if depth_structure == 1
        z0 = height./cos(deg2rad(gaze_angle));
        middle_dist = z0;
    end                                  %metersv
    cloud_dist                      = [middle_dist*.85, middle_dist*1.15];  % m range of distances   

    dot_density                     = .4;
    jitter                          = 2;                                    %deg to jitter by
    numSame                         = 5;                                   % number of apertures that have same texture
    cut                             = 4;                                     %cut off freq 1:4 (depends on speed, larger speed have lower cyc/deg corresponding to freq 1
    
    %-----Staircase Settings
    stairFrac                       = 10;                                    % fraction of distance between middle velcity and limit velocity
    depth_limits                    = [cloud_dist(1)/2, middle_dist];       %limits of depths to test, first is original deviation distance

    %-----Element Settings
    text_diam                       = 50;                                     % Arcmin 
    bps                             = text_diam+1;                          %so the moving texture is centered
    spatial_envelope                = 2;                                    % 0 = disk, 1 = Gabor, 2 = raised cosine
    devPos                          = [10; 6];                             % position of deviation (y-axis flipped when displaying)    

    exclude                         = [-2, -1, 2, 1];                          % Don't place elements where the FOE will be

    H_ecc_fix                       = 0;                                    % Horizontal fixation ecc (degs, neg is left)
    V_ecc_fix                       = 0;                                    % Vertical fixation ecc (degs, neg is up)

    experiment_id                   = 'pattern_detection';                  % Used in group filename
    fast                            = 1;                                    % Automatically trigger trials
    ITI                             = 1;                                    % Intertrial Inverval, Seconds
    background                      = 127;                                  % Grayscale Units
    fixate                          = 1;                                    % Present fixation spot during motion

    %-----Housekeeping----------------------
    % Scale things based on viewing distance, and convert other stuff to 
    % the units Psychtoolbox wants...
    tme                             = clock;
    logname = strcat(data_path,initials,'_',experiment_id, '_', date);

    text_pixelW                     = ceil(text_diam/scale_factor); %in pixels
    H_ecc_fix                       = H_ecc_fix*60/scale_factorX;
    V_ecc_fix                       = V_ecc_fix*60/scale_factorY;
    mv_length                       = ceil(stimulus_duration*frame_rate);

    %-----Spatial Envelope------------------
    if spatial_envelope == 1
        circle = 255*(1-ggaus(text_pixelW,.25*text_pixelW));
    elseif spatial_envelope == 2
        circle = 255*(1-coswin(text_pixelW,.3*text_pixelW,.5*text_pixelW));
    end
    
    %-----Open Screens----------------------
    Screen('Preference', 'SkipSyncTests', 1);
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    screens=Screen('Screens');
    screenNumber=max(screens);

    w=Screen('OpenWindow',screenNumber,background,[],[],2);
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    screen_rect = Screen('Rect',w);

    Screen('FillRect',w, background);
    Screen('Flip', w);
    Screen('FillRect',w, background);
    Screen('TextSize',w,20);

    if linearize
        calibration = load('0003_alta_141007.mat');
        table = calibration.calib.table;
        psychtoolbox_calib = repmat(table, 1, 3);
        Screen('LoadNormalizedGammaTable',screenNumber,psychtoolbox_calib);
    end 
    %-----Screen Landmarks------------
    sr_hor = round(screen_rect(3)/2); % Middle of the screen, horizontally, in pixels
    sr_ver = round(screen_rect(4)/2); % Middle of the screen, vertically, in pixels
    fix_hor = sr_hor+H_ecc_fix;     % Horizontal location of fixation cross, in pixels
    fix_ver = sr_ver+V_ecc_fix;     % Vertical location of fixation cross, in pixels
    movie_rect= [0,0,bps,bps];

    n_conditions = 0;
    for j=1:size(translation,1)
        for l = 1:length(theta)
            for k=1:n_staircases
                n_conditions = n_conditions+1;
                
                cond(n_conditions).rotation = theta(l);
                cond(n_conditions).LR = round(rand(1,n_trials)); %will the deviation be on the left or right, L=0, R=1
                cond(n_conditions).translate = translation(j,:);

                cond(n_conditions).count = 0; %what trial are we on
                cond(n_conditions).resp_history = zeros(1,n_trials); % whether you were right or wrong
                cond(n_conditions).n_correct = 0; % number of correct in a row - keeps track to change staircase
                cond(n_conditions).n_flips = 0;
                
                cond(n_conditions).step_idx = 1; % difficulty level (which velocity to choose)
                cond(n_conditions).step_history = zeros(1, n_trials); % to display staircase
                cond(n_conditions).steps = zeros(2,stairFrac); % list of possible velocities for the staircase
                cond(n_conditions).dev_history = zeros(2,n_trials); % what deviation vector was it
                cond(n_conditions).dev_angle_history = zeros(1,n_trials); % history of angles
                cond(n_conditions).dev_speed_history = zeros(1,n_trials); % percent speed
       
                cond(n_conditions).velocity_range = zeros(2,2); %range of possible velocities base on range of depths
                cond(n_conditions).angle_btwn_constraints = theta; %calculate dot product between on constraint and off constraint
                cond(n_conditions).z0_middle_dist = [0 0]; %save fixation depth and middle distance of cloud

            end
        end
    end
    results = zeros(n_conditions,n_trials);
    
    % calculate all possible velocities for the deviations and save in
    % condition struct
    velocity_steps = zeros(2*size(translation,1), stairFrac);
    velocity_steps_off = zeros(size(velocity_steps));
    velocity_range = zeros(size(translation,1), 2);
    velocity_range_off = zeros(2*size(translation,1), 2);
    angle_btwn_constraints = zeros(2,size(translation,1));
    
    % get velocity steps vv
    for c = 1:length(cond)


        end_vel = calculate_cloud_flow(depth_limits(2), [devPos(1); devPos(2)], cond(c).translate, view_dist, z0);
        %         end_vel = [0; -2.5125];
        start_vel = deviation_scale.* end_vel;

        % get steps
        segment = end_vel-start_vel; %vector describing segment between starting velocity, and middle velocity
        steps = [linspace(0, segment(1), stairFrac); linspace(0, segment(2), stairFrac)]; %amount to move along segment
        velocity_steps = start_vel+steps;
        
         % rotate by theta
        velocity_steps = get_theta_off_constraint(velocity_steps, end_vel, cond(c).rotation);  
        
        % assign velocities to appropriate condition
        cond(c).steps = velocity_steps;
        
        % velocity range
        velocity_range(:,1) = calculate_cloud_flow(cloud_dist(1), devPos, cond(c).translate, view_dist, z0);
        velocity_range(:,2) = calculate_cloud_flow(cloud_dist(2), devPos, cond(c).translate, view_dist, z0);
        cond(c).velocity_range = velocity_range;
    end
    
    %-----Randomize Trials------------
    total_trials = n_trials*n_conditions;
    perm = randperm(total_trials);
    perm = mod(perm,n_conditions)+1;
    duration_check = zeros(size(perm));

    Screen('DrawText',w,'Heading discrimination: Pattern',500,240,250);
    Screen('DrawText',w,'Use Left/Right Arrows to respond',500,270,250);
    Screen('DrawText',w,'press V to start',500,300,250);
    Screen('Flip',w);

    if mouseclick
        validKey = 0;
        while ~validKey
            [clicks,x,y,whichButton] = GetClicks();
            find(clicks)
            if clicks >= 1
                validKey = 1;
            end
        end
    else
        FlushEvents('keyDown');
        validKey = 0;
        while ~validKey
            [secs, keyCode, deltaSecs] = KbWait();
            find(keyCode)
            if keyCode(KbName('v'))
                validKey = 1;
            end
        end
    end

    Screen('FillRect', w, background);
    Screen('Flip', w);
    tic;

    if test
        total_trials = 1;
    else
    end

    for trial = 1:total_trials
        aa = GetSecs; 

        % Draw the white fixation cross
        Screen('FillRect',w, background);
        Screen('DrawDots', w, [fix_hor; fix_ver], 4, 250, [], 2); 
        Screen('Flip',w);
        
        if test
            perm(trial) = 1;
        end

        cond(perm(trial)).count = cond(perm(trial)).count+1;
        % Set up the world

%         translate = cond(perm(trial)).translate;
%         if depth_structure == 1
%             [dots_m, dots_deg, dev0, dev1] = make_dot_plane(dot_density, jitter, bumpiness,height, gaze_angle, wall, view_window, text_diam, exclude, devPos);
%         elseif depth_structure == 2
%         else
%             [dots_m, dots_deg, dev0, dev1] = make_dot_cloud(dot_density, jitter, cloud_dist, view_window, text_diam, exclude, devPos);
%         end

        % repeat same velocity in surround
        [dots_deg, dev0, dev1] = make_surround(dot_density, jitter, view_window, devPos); 
        velocity_field = repmat(cond(1).steps(:,end), 1, length(dots_deg));
        
        
        if scramble
            translate = cond(perm(trial)).translate;
%           % scramble optic flow info
            [dots_m, dots_deg, dev0, dev1] = make_surround_withflow(dot_density, jitter, cloud_dist, view_window, devPos);
            velocity_field = calculate_cloud_flow(dots_m(3,:), dots_deg, translate, view_dist, z0);

            %scramble left and right
            
            shortened = velocity_field;
            
            shortened(:,dev0) = [];
            shortened(:,dev1) = [];
            leftidx = find(shortened(1,:) <0);
            rightidx = find(shortened(1,:) >0);
            scrambledidxleft = randsample(leftidx, length(shortened)/2);
            scrambledidxright = randsample(rightidx, length(shortened)/2);
            scrambled = shortened;
            scrambled(:,leftidx) = scrambled(:, scrambledidxleft);
            scrambled(:,rightidx) = scrambled(:, scrambledidxright);
        
            % put centers back
            full_scramble = [scrambled(:,1:dev0-1) velocity_field(:,dev0) scrambled(:,dev0:dev1-2) velocity_field(:, dev1) scrambled(:, dev1-1:end)];
        
            %save as velocity field
            velocity_field = full_scramble;
        end
        
        if cond(perm(trial)).LR(cond(perm(trial)).count) %if it's on the right
            idx = dev1;
            velocity_field(:,idx) = cond(perm(trial)).steps(:,cond(perm(trial)).step_idx);
        else % if it's on the left
            idx = dev0;
            velocity_field(:,idx) = cond(perm(trial)).steps(:,cond(perm(trial)).step_idx);
            if scramble
                velocity_field(:,idx) = [-cond(perm(trial)).steps(1,cond(perm(trial)).step_idx); cond(perm(trial)).steps(2,cond(perm(trial)).step_idx)];
            end
        end
        

        %pixels/s
        velocity_field_pix = velocity_field*60*1/scale_factor;
        %convert to pixels/frame
        velocity_field_ppf = velocity_field_pix/frame_rate;
        
        %-----1/f texture-----------------------
        sz = max(vecnorm(velocity_field_pix))*2*stimulus_duration + text_pixelW;

        numDrifts = ceil(length(velocity_field)/numSame);
        tex = zeros(1,numDrifts);
        for d = 1:numDrifts
            drift = oneoverfcut(1, ceil(sz), ceil(sz), cut);
            tex(d)=Screen('MakeTexture', w, 255*drift);
        end
        drift_center = [size(drift,1)/2, size(drift,1)/2];
        which_tex = randi(numDrifts, 1,length(velocity_field)); %assign one of drifts to each patch for the movie

        windowTex = Screen('MakeTexture',w,cat(3,128*ones(text_pixelW),255*(1-coswin(text_pixelW,.3*text_pixelW,.5*text_pixelW))));
        
        % get positions to draw rectangles
        pos0 = CenterRectOnPoint(movie_rect, drift_center(1), drift_center(1));
        pos = CenterRectOnPoint(movie_rect,dots_deg(1,:)'*60/scale_factorX+sr_hor,dots_deg(2,:)'*60/scale_factorY+sr_ver)'; %big Y numbers are lower on the screen in pixels so we have to flip

        source = cell(1,length(stimulus_duration*frame_rate));
        for tt = 1:ceil(stimulus_duration*frame_rate)
            for s = 1:length(velocity_field)
                source{tt}(:,s) = [pos0(1)-(tt-1)*velocity_field_ppf(1,s); pos0(2)-(tt-1)*velocity_field_ppf(2,s); ...
                    pos0(3)-(tt-1)*velocity_field_ppf(1,s); pos0(4)-(tt-1)*velocity_field_ppf(2,s)];
            end
        end
        
        %Finish the ITI
        WaitSecs(ITI-(GetSecs-aa));

        % Draw the red fixation cross 
        Screen('DrawDots', w, [fix_hor; fix_ver], 4, [250 0 0], [], 2);
        Screen('Flip',w);

    %     FlushEvents('keyDown'); 
    %     priorityLevel=MaxPriority(w); 
    %     Priority(priorityLevel);

        % Play the movie
        
        StimulusOnsetTime = zeros(1,mv_length);
        aa = GetSecs;
        tag = clock;

        for frame = 1:mv_length 
            Screen('DrawTextures', w, tex(which_tex),source{frame},pos, 0);
            Screen('DrawTextures', w, windowTex, [], pos, 0); 
            if fixate 
                Screen('DrawDots', w, [fix_hor; fix_ver], 4, [250 0 0], [], 2);
            end
%             [VBLTimestamp, StimulusOnsetTime(frame), FlipTimestamp,
%             Missed, Beampos] = Screen('Flip',w);
            Screen('Flip', w);
        end

        duration_check(trial) = GetSecs-aa;
        Screen('FillRect',w, background);
        if fixate
            Screen('DrawDots', w, [fix_hor; fix_ver], 4, [250 0 0], [], 2);
        end
        Screen('Flip',w); 

        if mouseclick
            validKey = 0;
            while ~validKey
                [clicks,x,y,whichButton] = GetClicks([],[0],[]);
        %         key = KbName();
        %         if key(1) == 'E'
        %             Screen('CloseAll');
        %             Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
        %             Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
        %             ListenChar(1);
        %             ShowCursor;
        %             break
                if whichButton == 1 
                    validKey = 1;
                    resp = 0;
                elseif whichButton == 2
                    validKey = 1;
                    resp = 1;
                else
                    Beeper('low');
                end
            % Provides feedback during training phase
                if trainingMode
                    if cond(perm(trial)).LR(cond(perm(trial)).count) == resp 
                        Screen('DrawText',w,'Correct',625,375,250);
                        Screen('Flip',w)
                        WaitSecs(ITI)
                    else
                        Screen('DrawText',w,'Incorrect',625,375 ,250);
                        Screen('Flip',w)
                        WaitSecs(ITI)  
                    end
                end
            end
        else
     %     Get the response
            validKey = 0;
            while ~validKey
                [secs, keyCode, deltaSecs] = KbWait();
                if keyCode(KbName('ESCAPE'))
                    Screen('CloseAll');
                    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
                    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
                    ListenChar(1);
                    ShowCursor;
                    break
                elseif keyCode(KbName('LeftArrow'))
                    validKey = 1;
                    resp = 0;
                elseif keyCode(KbName('RightArrow')) 
                    validKey = 1;
                    resp = 1;
                else
                    Beeper('low');
                end
%            Provides feedback during training phase
                if trainingMode
                    if cond(perm(trial)).LR(cond(perm(trial)).count) == resp 
                        Screen('DrawText',w,'Correct',625,375,250);
                        Screen('Flip',w)
                        WaitSecs(ITI)
                    else
                        Screen('DrawText',w,'Incorrect',625,375 ,250);
                        Screen('Flip',w)
                        WaitSecs(ITI)  
                    end
                end
            end
        end
    %     Priority(0); 
    %     % Log the trial
    %     if ~trainingMode
    %         fid = fopen(logname,'a');
    %         fprintf(fid,'%s, %s, %d, %d, %d, %d, %d, %f, %d, %f, %f, %f , %d, %d, %d, %f,%f, %f %f\n', ...
    %                 initials,experiment_id, tag(1), tag(2), tag(3), tag(4), tag(5), tag(6), ...
    %                 trial, cond(perm(trial)).translate(1), cond(perm(trial)).translate(2), cond(perm(trial)).translate(3), ...
    %                 cond(perm(trial)).dev_index, cond(perm(trial)).LR(cond(perm(trial)).count), resp, percentSpeed, diffAngle, ...
    %                 stimulus_duration, duration_check(trial));        
    %         fclose(fid); 
    %     end

    % Tell the staircase what happened
    %report just rhs staircase
    cond(perm(trial)).dev_history(:,cond(perm(trial)).count) = cond(perm(trial)).steps(:,cond(perm(trial)).step_idx);
    cond(perm(trial)).step_history(cond(perm(trial)).count) = cond(perm(trial)).step_idx;
    
    % percent speed
    cond(perm(trial)).percentSpeed(cond(perm(trial)).count) = norm(cond(perm(trial)).steps(:,cond(perm(trial)).step_idx))/norm(cond(:,perm(trial)).steps(:,end));
    
    % angle difference
    unitDev = cond(perm(trial)).steps(:,cond(perm(trial)).step_idx)/norm(cond(perm(trial)).steps(:,cond(perm(trial)).step_idx));
    unitMiddle = cond(perm(trial)).steps(:,end)/norm(cond(perm(trial)).steps(:,end));
    cond(perm(trial)).dev_angle_history(cond(perm(trial)).count) = abs(rad2deg(acos(unitDev'*unitMiddle)));
    
    % if they haven't flipped yet
    if cond(perm(trial)).n_flips < 1
        % if they're correct
        if cond(perm(trial)).LR(cond(perm(trial)).count) == resp
            cond(perm(trial)).step_idx = cond(perm(trial)).step_idx + 2; % make it harder by two steps
            cond(perm(trial)).resp_history(cond(perm(trial)).count) = 1; % report a correct response in history
        else %if they're wrong
            cond(perm(trial)).step_idx = cond(perm(trial)).step_idx - 1; % maket it easier by 1 step
        end
        if cond(perm(trial)).count > 1 && cond(perm(trial)).resp_history(cond(perm(trial)).count) ~= cond(perm(trial)).resp_history(cond(perm(trial)).count-1)
            %if after trial one and the last answer is different than the
            %current on, increase n_flip
            cond(perm(trial)).n_flips = 1;
        end
    else %if there has been a flip
        if cond(perm(trial)).LR(cond(perm(trial)).count) == resp %if they are correct
            cond(perm(trial)).n_correct = cond(perm(trial)).n_correct+1;
            cond(perm(trial)).resp_history(cond(perm(trial)).count) = 1; % report a correct response in history
            if cond(perm(trial)).n_correct == 2
                cond(perm(trial)).step_idx = cond(perm(trial)).step_idx + 1; % make it harder by a step
                cond(perm(trial)).n_correct = 0; % reset number correct
            else
                %dont change step idx, keep counting how many correct in a row
            end
        else %if they are wrong
            cond(perm(trial)).n_correct = 0; % reset number correct in a row
            cond(perm(trial)).step_idx = cond(perm(trial)).step_idx - 1; % make it easier by a step
        end
    end
    
    % check: if the staircase goes out of bounds/preset velocities reset to
    % nearest index
    if cond(perm(trial)).step_idx > stairFrac
        cond(perm(trial)).step_idx = stairFrac;
    elseif cond(perm(trial)).step_idx < 1
        cond(perm(trial)).step_idx = 1;
    end

    end
    clc;

    sca;
    today = date;
    day = today(1:2);
    month = today(4:6);
    year = today(end-3: end);
    if ~trainingMode
        fname = sprintf ( '%s-%s-%s_%s_%s.mat', year, month, day,initials, session);
        folder = [data_path fname];
        save(folder, 'cond')
    end
% catch ME
%     %this "catch" section executes in case of an error in the "try" section
%     %above.  Importantly, it closes the onscreen window if its open.
%     ListenChar(1);
%     ShowCursor;
%     Screen('CloseAll');
%     Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
%     Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
%     Priority(0);v
%v end %try..catch   