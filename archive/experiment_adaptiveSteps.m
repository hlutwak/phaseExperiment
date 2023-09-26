
%% Phase experiment, but staircasing with linear changes along the constraint
% updated HL 9/6/2020
clear all;close all;clc;
addpath(genpath('/Applications/Psychtoolbox'))
%vvvvv 
    %-----Subject Settings------------------ 
    initials                        = 'HL';                                 % Subject initials used to save file
    session                         = 00';
    trainingMode                    = 1;                                    % 1 = gives feedback after each trial (training), 0 = no feedback (actual experiment)

    %-----Rig Settings----------------------
    mouseclick                      = 1;                                    % use mouseclick instead of keys
    view_dist                       = .35;                                  % m .57; psychophysics room: .35 .50
    screensize                      = [.405 .298];                           % m
    pixels                          = [1600 1200];
    frame_rate                      = 85;   

    view_window                     = round([rad2deg(atan(screensize(1)/2/view_dist)) rad2deg(atan(screensize(2)/2/view_dist))]);  % X,Y centered around fixation [60 46] [54 40];                                                                      % [36 27] for laptop, [55 32] for monitor
    scale_factor                    = sqrt((view_window(1)*60*view_window(2)*60)/(pixels(1)*pixels(2)))*1;          % Arcmin/pixel 1.78 2.25 1.92dell: 1680x1050, psychophysics room: 1600x1200                                % Screen frame rate (hz) psychophysics room: 85
    scale_factorX = view_window(1)*60/pixels(1);
    scale_factorY = view_window(2)*60/pixels(2);
    linearize                       = 0;                                    % Use calibrated LUT (do this when available)                    

    %-----Trial Settings------------------ 
    subjectNumber                   = '01';                                 % Number used to save eyelink file
    n_trials                        = 20;                                  % Number of trials
    translation                     = [ 0, .03, 1];              
%                                       [ 0, -0.15, .65]];                 % deg/s *up is down and down is up
    depth_structure                 = 1;                                    % 0 = random cloud, 1 = ground plane, 2 = hallway  
    if depth_structure
        gaze_angle                  = 77;                                   %deg from staring straight at ground
        translation                 = getT_eyecoords(90-gaze_angle, translation);
        height                      = 1.65;                                 % height in meters
        bumpiness                   = 0;                                  % random noise in ground plane
        wall                        = 0;                                    % ground plane + wall starting in upper visual field
    end
    rotate_constraint               = 1;                                    % get constraint line by rotating by theta rather than mirroring
    if rotate_constraint
        theta                       = 270;
    end
    fixateMovement                  = 0;                                    % 1 = red fixation dot moves, 0 = stationary red fixation dot
    eyeTracking                     = 0;                                    % 1 = eye tracking on, 0 = eye tracking off
    linearize                       = 0;                                    % Use calibrated LUT (do this when available)
    data_path                       = './data/';
    stimulus_duration               = .5;                                    % seconds
    
    test                            = 1;                                    % just shows a few stimv
    if test
        stimulus_duration = 2;
    end
    n_staircases                    = 2;                                    % Want to run multiple staircases

    % % Decides if experiment is simulated rotation or eye tracking and sets variables accordingly
    % linear_velocity                 = (rotation * .0174533) * 30;           % m/s, calculated from rotation (angular velocity)
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
    end                                  %meters
    cloud_dist                      = [middle_dist*.85, middle_dist*1.15];  % m range of distances   

    dot_density                     = .4;
    jitter                          = 2;                                    %deg to jitter by
    %-----Staircase Settings
    stairFrac                       = 10;                                    % fraction of distance between middle velcity and limit velocity
    depth_limits                    = [cloud_dist(1)/2, middle_dist];       %limits of depths to test, first is original deviation distance

    %-----Element Settings
    gabor_diam                      = 66;                                  % Arcmin 
    gabor_sf                        = 2.2;                                    % c/deg   
    gabor_contrast                  = 100;                                  % percent, for the pattern; each element will be half this
    spatial_envelope                = 2;                                    % 0 = disk, 1 = Gabor, 2 = raised cosine
    devPos                          = [10; 7];                             % position of deviation (y-axis flipped when displaying)    

    exclude                         = [-5 -3 5 3];                          % Don't place elements where the FOE will be

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

    gabor_size                      = gabor_diam/scale_factor;
    stimulus_radius                 = round(gabor_size/2);
    H_ecc_fix                       = H_ecc_fix*60/scale_factorX;
    V_ecc_fix                       = V_ecc_fix*60/scale_factorY;
    mv_length                       = ceil(stimulus_duration*frame_rate);
    f                               = (gabor_sf*scale_factor/60)*2*pi;
    angle                           = 0;
    a                               = cos(angle)*f;
    b                               = sin(angle)*f;
    amplitude                       = background;

    %-----Spatial Envelope------------------
    [x,y]=meshgrid(-stimulus_radius:stimulus_radius,-stimulus_radius:stimulus_radius);
    bps = (stimulus_radius)*2+1;
    circle=((stimulus_radius)^2-(x.^2+y.^2));
    for ii=1:bps; for j =1:bps; if circle(ii,j) < 0; circle(ii,j) = 0; else circle(ii,j) = 1; end; end;
    end;
    if spatial_envelope == 1
        circle = (exp(-(((x)/(sqrt(2)*Gaussian_stdev/6)).^2)-((y/(sqrt(2)*Gaussian_stdev/2)).^2)).*circle);
    elseif spatial_envelope == 2
        R = (sqrt(x.^2 + y.^2) + eps).*circle;
        R = R/max(max(R));
        cos2D = (cos(R*pi)+1)/2;
        circle = (cos2D.*circle);
    end
    circle = circle*255/2;

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

    %-----Make the grating in all phases-----
    gabor = zeros(1,360);
    for ii=1:360
        grating = round(((sin(a*x+b*y+ii*pi/180)*amplitude)+background));
        gabor(ii) = Screen('MakeTexture',w,cat(3,grating,circle));
    end

    n_conditions = 0;
    for j=1:size(translation,1)
        for l = 1:2
            for k=1:n_staircases
                n_conditions = n_conditions+1;
                
                cond(n_conditions).dev_index = mod(l,2); %dev_index = 1 for off constraint, 0 for on constraint
                cond(n_conditions).LR = round(rand(1,n_trials)); %will the deviation be on the left or right, L=0, R=1
                cond(n_conditions).translate = translation(j,:);
                cond(n_conditions).contrast = gabor_contrast/100;
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
                cond(n_conditions).angle_btwn_constraints = 0; %calculate dot product between on constraint and off constraint
                cond(n_conditions).z0_middle_dist = [0 0]; %save fixation depth and middle distance of cloud

            end
        end
    end
    results = zeros(n_conditions,n_trials);
    
    % calculate all possible velocities for the deviations and save in
    % condition struct
    velocity_steps = zeros(2*size(translation,1), stairFrac);
    velocity_steps_off = zeros(size(velocity_steps));
    velocity_range = zeros(2*size(translation,1), 2);
    velocity_range_off = zeros(2*size(translation,1), 2);
    angle_btwn_constraints = zeros(2,size(translation,1));
    
    % get velocities for on constraint
    for ii = 1:size(translation,1)
        if depth_structure == 1
            z_dev = raytrace2eyeZ(devPos, height, gaze_angle);
            depth_limits                    = [z_dev-3*bumpiness, z_dev];       %limits of depths to test, first is original deviation distance
            
            cloud_dist(1) = z_dev-1*bumpiness; %range of velocities possible at this location based on noisy bumpy ground
            cloud_dist(2) = z_dev+1*bumpiness;
        end
        start_vel = calculate_cloud_flow(depth_limits(1), devPos, translation(ii,:), view_dist, z0);
        end_vel = calculate_cloud_flow(depth_limits(2), devPos, translation(ii,:), view_dist, z0);

        velocity_range(2*ii-1:2*ii,1) = calculate_cloud_flow(cloud_dist(1), devPos, translation(ii,:), view_dist, z0);
        velocity_range(2*ii-1:2*ii,2) = calculate_cloud_flow(cloud_dist(2), devPos, translation(ii,:), view_dist, z0);
        
        % get steps
        segment = end_vel-start_vel; %vector describing segment between starting velocity, and middle velocity
        steps = [linspace(0, segment(1), stairFrac); linspace(0, segment(2), stairFrac)]; %amount to move along segment
        velocity_steps(2*ii-1:2*ii,:) = start_vel+steps;
        
        if rotate_constraint
             % calculate for off the constraint by rotating by theta
            velocity_steps_off(2*ii-1:2*ii, :) = get_theta_off_constraint(velocity_steps(2*ii-1:2*ii,:), end_vel, theta);
            velocity_range_off(2*ii-1:2*ii, :) = get_theta_off_constraint(velocity_range(2*ii-1:2*ii,:), end_vel, theta);
    
        else
            % calculate for off the constraint
            velocity_steps_off(2*ii-1:2*ii, :) = get_off_constraint(velocity_steps(2*ii-1:2*ii,:), end_vel);
            velocity_range_off(2*ii-1:2*ii, :) = get_off_constraint(velocity_range(2*ii-1:2*ii,:), end_vel);
        end
        
        
        % calculate dot product between on/off constraint lines
        onDir = velocity_range(2*ii-1:2*ii, 1) - velocity_range(2*ii-1:2*ii, 2);
        normOn = onDir/norm(onDir);
        offDir = velocity_range_off(2*ii-1:2*ii, 1) - velocity_range_off(2*ii-1:2*ii, 2);
        normOff = offDir/norm(offDir);
        angle_btwn_constraints(ii) = rad2deg(acos(normOn'*normOff));
    end
    
    % assign appropriate velocities to staircases
    
    for ii = 1: length(cond)
        if mean(cond(ii).translate == translation(1,:))== 1 %if first translation
            if cond(ii).dev_index % 1 = on constraint
                cond(ii).steps = velocity_steps(1:2, :);
                cond(ii).velocity_range = velocity_range(1:2,:);
            else
                cond(ii).steps = velocity_steps_off(1:2,:);
                cond(ii).velocity_range = velocity_range_off(1:2,:);
            end
            cond(ii).angle_btwn_constraints = angle_btwn_constraints(1);
        else
            if cond(ii).dev_index
                cond(ii).steps = velocity_steps(3:4,:);
                cond(ii).velocity_range = velocity_range(3:4,:);
            else
                cond(ii).steps = velocity_steps_off(3:4,:);
                cond(ii).velocity_range = velocity_range_off(3:4,:);
            end
            cond(ii).angle_btwn_constraints = angle_btwn_constraints(2);
        end
        cond(ii).z0_middle_dist = [z0, middle_dist];
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

        cond(perm(trial)).count = cond(perm(trial)).count+1;
        % Set up the world
        angle_all = [];
        speed_all = [];
        pos_all = [];

        translate = cond(perm(trial)).translate;
        if depth_structure == 1
            [dots_m, dots_deg, dev0, dev1] = make_dot_plane(dot_density, jitter, bumpiness,height, gaze_angle, wall, view_window, gabor_diam, exclude, devPos);
        elseif depth_structure == 2
        else
            [dots_m, dots_deg, dev0, dev1] = make_dot_cloud(dot_density, jitter, cloud_dist, view_window, gabor_diam, exclude, devPos);
        end
        
        velocity_field = calculate_cloud_flow(dots_m(3,:), dots_deg, translate, view_dist, z0);

        if cond(perm(trial)).LR(cond(perm(trial)).count) %if it's on the right
            idx = dev1;
            velocity_field(:,idx) = cond(perm(trial)).steps(:,cond(perm(trial)).step_idx);
        else % if it's on the left
            idx = dev0;
            velocity_field(:,idx) = [-cond(perm(trial)).steps(1,cond(perm(trial)).step_idx); cond(perm(trial)).steps(2,cond(perm(trial)).step_idx)];
            
        end
        
        nGabors(ii) = size(dots_deg,2);  
        posA = CenterRectOnPoint(movie_rect,dots_deg(1,:)'*60/scale_factorX+sr_hor,dots_deg(2,:)'*60/scale_factorY+sr_ver)'; %big Y numbers are lower on the screen in pixels so we have to flip
        pos = [posA posA];    

        speed_ver = velocity_field(2,:);    
        speed_hor = velocity_field(1,:);

        angle_ver = 270*(speed_ver>=0)+90*(speed_ver<0); % in deg  
        angle_hor = 180*(speed_hor>0); % in deg
        speed_ver = abs(velocity_field(2,:));
        speed_hor = abs(velocity_field(1,:));

        speed = zeros(1,length(speed_ver)*2); 
        angle = speed;
        speed = [speed_hor speed_ver];
        angle = [angle_hor angle_ver];
        angle_all = [angle_all angle];
        speed_all = [speed_all speed];
        pos_all = [pos_all pos];

        angle = angle_all;   
        speed = speed_all; 
        pos = pos_all;
        phase_step = round(speed*gabor_sf*360/frame_rate)';

        phases = repmat(ceil(rand(size(speed,2),1)*360),1,mv_length)+repmat(phase_step,1,mv_length).*repmat(0:(mv_length-1),size(speed,2),1);
        phases = rem(phases,360);
        phases = phases + (phases==0)*360; 

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
            Screen('DrawTextures', w, gabor(phases(:,frame)),movie_rect,pos,angle); 
            if fixate 
                Screen('DrawDots', w, [fix_hor; fix_ver], 4, [250 0 0], [], 2);
            end
            [VBLTimestamp, StimulusOnsetTime(frame), FlipTimestamp, Missed, Beampos] = Screen('Flip',w);
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
            if cond(perm(trial)).n_correct == 2
                cond(perm(trial)).step_idx = cond(perm(trial)).step_idx + 1; % make it harder by a step
                cond(perm(trial)).resp_history(cond(perm(trial)).count) = 1; % report a correct response in history
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
    fname = sprintf ( '%s-%s-%s_%s_%s.mat', year, month, day,initials, session);
    folder = [data_path fname];
    save(folder, 'cond')
    
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

%% 

%         
% figure(1)
% hold on
% scatter(dots_deg(1,:), -dots_deg(2,:), 'b')   
% hold on, scatter(dots_deg(1,idx), -dots_deg(2,idx), 'r') 
% axis equal
% 
% figure(2) %plot deviations
% % plot velocity at middle depth
% hold on, quiver(0,0,cond(perm(trial)).original(1), -cond(perm(trial)).original(2),'LineWidth', 3, 'color', 'k', 'AutoScaleFactor',1, 'DisplayName', 'middle depth')
% % plot deviation off constraint
% hold on, quiver(0,0,deviation(1), -deviation(2), '--','LineWidth', 3,'color', [0.8500, 0.3250, 0.0980], 'AutoScaleFactor',1, 'DisplayName', 'off constraint')
% 
% hold on, quiver(0,0,velocity_field_dev(1), -velocity_field_dev(2),'--','LineWidth', 3, 'color', 'b', 'AutoScaleFactor',1, 'DisplayName', 'on constraint')
% hold on, quiver(0,0,velocity_field(1,idx), -velocity_field(2,idx), 'color','r', 'AutoScaleFactor',1, 'DisplayName', 'chosen deviation')
% 
% %plot range of velocities
% hold on, plot(cond(1).velocity_range(1,:), -cond(1).velocity_range(2,:),'color',[.75 0 0], 'DisplayName', 'off constraint range')
% hold on, plot(cond(2).velocity_range(1,:), -cond(2).velocity_range(2,:), 'color',[0 0 .75], 'DisplayName', 'on constraint range')
% sca
% 
% axis equal 
% xlim([-10,10]), ylim([-10,10]) 
% legend  