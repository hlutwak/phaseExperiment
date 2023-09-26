%% Figures
figure
scatter(dots_deg(1,:), -dots_deg(2,:), 'b')   
hold on, scatter(dots_deg(1,idx), -dots_deg(2,idx), 'r') 
axis equal

figure, quiver(0,0,original(1), -original(2),'LineWidth', 2, 'color', 'k')
% hold on, quiver(0,0,deviation(1), -deviation(2), '--','LineWidth', 2,'color', 'c')
hold on, quiver(0,0,velocity_field_dev(1), -velocity_field_dev(2),'--','LineWidth', 2, 'color', 'b')
axis equal, xlim([-5,5]), ylim([-5,5]) 
legend('original velocity', 'off constraint', 'along constraint', 'chosen deviation') 

%%
figure
scatter(dots_deg(1,:), -dots_deg(2,:), 'b')   
hold on, scatter(dots_deg(1,idx), -dots_deg(2,idx), 'r') 
axis equal
figure, quiver(0,0,original(1), -original(2),'LineWidth', 3, 'color', 'k')
hold on, quiver(0,0,deviation(1), -deviation(2), 'LineWidth', 3,'color', c(3,:))
hold on, quiver(0,0,velocity_field_dev(1), -velocity_field_dev(2),'LineWidth', 3, 'color', c(7,:))
 
hold on, quiver(0,0,velocity_field(1,idx), -velocity_field(2,idx), 'color','r')
axis equal, xlim([-5,5]), ylim([-5,5]) 
legend('original velocity', 'off constraint', 'along constraint', 'chosen deviation')   

%% line segment plot
% clear all;clc;
addpath(genpath('/Applications/Psychtoolbox'))
 
%-----Settings to Modify------------------ 
initials                        = 'test';                          % Subject initials used to save file
session                         = 17;
subjectNumber                   = '00';                                 % Number used to save eyelink file
n_trials                        = 25;                                   % Number of trials  per staircase
translation                     = [ 0, -0.01, 1.9021];                 
%                                    [ 0, -0.01, 1.9021]];              % deg/s *up is down and down is up
fixateMovement                  = 0;                                    % 1 = red fixation dot moves, 0 = stationary red fixation dot
eyeTracking                     = 0;                                    % 1 = eye tracking on, 0 = eye tracking off
test                            = 0;                                    % just shows a few stimv
trainingMode                    = 0;                                    % 1 = gives feedback after each trial (training), 0 = no feedback (actual experiment)
linearize                       = 0;                                    % Use calibrated LUT (do this when available)
% data_path                       ='/ser/1.1/p1/bonnen/Documents/MATLAB/phaseExperiment/phaseExperiment/data/';         % Folder for saving data files
data_path                       = './data';
stimulus_duration               = 1;                                    % seconds
if test
    stimulus_duration = 10;
end
n_staircases                    = 1;                                    % Want to run multiple staircases
 
% % Decides if experiment is simulated rotation or eye tracking and sets variables accordingly
% linear_velocity                 = (rotation * .0174533) * 30;           % m/s, calculated from rotation (angular velocity)
% if fixateMovement == 1
%     rotation = 0;vv
% elseif fixateMovement == 0
%     linear_velocity = 0;
% end
 
%-----Experiment Settings, Don't Change----
%-----Array Settings
middle_dist                     = 20;                                   %meters
cloud_dist                      = [middle_dist*.75, middle_dist*1.25];  % m range of distances     
dot_density                     = .06;
 
%-----Staircase Settings
step                            = .5;                                    %meters
limits                          = [cloud_dist(1)*(1/2.5) middle_dist];
 
%-----Element Settings
gabor_diam                      = 75;                                  % Arcmin 
gabor_sf                        = 2;                                    % c/deg   
gabor_contrast                  = 100;                                  % percent, for the pattern; each element will be half this
spatial_envelope                = 2;                                    % 0 = disk, 1 = Gabor, 2 = raised cosine
devPos                          = [-12;-7];                             % position of deviation (y-axis flipped when displaying)    
 
view_window                     = [60 46];                              % X,Y centered around fixation [60 46];
                                                                        % [36 27] for laptop, [55 32] for monitor
exclude                         = [-8 -6 8 6];                       % Don't place elements where the FOE will be
 
 
H_ecc_fix                       = 0;                                    % Horizontal fixation ecc (degs, neg is left)
V_ecc_fix                       = 0;                                    % Vertical fixation ecc (degs, neg is up)
 
experiment_id                   = 'pattern_detection';                    % Used in group filename
fast                            = 1;                                    % Automatically trigger trials
ITI                             = 1;                                    % Intertrial Inverval, Seconds
background                      = 127;                                  % Grayscale Units
fixate                          = 1;                                    % Present fixation spot during motion
 
%-----Rig Settings----------------------
view_dist                       = .35;                                    % m .57; psychophysics room: .35
scale_factor                    = 2.25;                                 % Arcmin/pixel 1.78 2.25 dell: 1680x1050, psychophysics room: 1600x1200
frame_rate                      = 85;                                   % Screen frame rate (hz) psychophysics room: 85
linearize                       = 0;                                    % Use calibrated LUT (do this when available)
 
%-----Housekeeping----------------------
% Scale things based on viewing distance, and convert other stuff to 
% the units Psychtoolbox wants...
tme                             = clock;
logname = strcat(data_path,initials,'_',experiment_id, '_', date);
 
gabor_size                      = gabor_diam/scale_factor;
stimulus_radius                 = round(gabor_size/2);
H_ecc_fix                       = H_ecc_fix*60/scale_factor;
V_ecc_fix                       = V_ecc_fix*60/scale_factor;
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
 
%-----Set Up Conditions-----------
n_conditions = 0;
for j=1:size(translation,1)
        for k=1:n_staircases
            n_conditions = n_conditions+1; 
            cond(n_conditions).LR = ones(1,n_trials); %will the deviation be on the left or right, L=0, R=1
            cond(n_conditions).translate = translation(j,:); 
 
            cond(n_conditions).contrast = gabor_contrast/100;
 
            cond(n_conditions).deviation_depth = [15 25]; % visualizing different depths
        end
end

%-----Randomize Trials------------
total_trials = n_trials*n_conditions;
perm = randperm(total_trials);
perm = mod(perm,n_conditions)+1;
duration_check = zeros(size(perm));
trial = 1;
translate = translation;
figure
on = winter(length(cond.deviation_depth));
off = autumn(length(cond.deviation_depth));
difference = zeros(2,length(cond.deviation_depth));

diffAngle = zeros(length(cond.deviation_depth),1);
for d = 1: length(cond.deviation_depth)
    [dots_m, dots_deg] = make_dot_cloud(dot_density, cloud_dist, view_window, gabor_diam, exclude, devPos);
    velocity_field = calculate_cloud_flow(dots_m, dots_deg, translate, view_dist);
 
    idx = 2;
    dev_m = [dots_m(1:2,2); cond(perm(trial)).deviation_depth(d)]; % get deviation on right

    original = velocity_field(:,idx); 
    cond(perm(trial)).dev_original = original;
    %--get deviation ON CONSTRAINT---
    dev_deg =  dots_deg(:,idx);
    velocity_field_dev = calculate_cloud_flow(dev_m, dev_deg, translate, view_dist);
 
    %-- get deviation same change in speed and  direction OFF CONSTRAINT-- 
    unitOriginal = velocity_field(:,idx)/norm(velocity_field(:,idx));
    theta = atan(unitOriginal(2)/unitOriginal(1)); %get angle between unit original velocity and x-axis
    R = [cos(theta) sin(theta); -sin(theta) cos(theta)]; % rotates based on tan(theta) 
    rotated = R*velocity_field_dev; %rotate the deviated velocity along constraint
    rotated(2) = -rotated(2); %reflect over x-axis (rotated direction of original velocity)
    Rback = R'; %rotation matrix to go back  
    deviation = Rback*rotated;   
 
    %-- calculate the difference
    difference(:,d) = velocity_field_dev - original;
    percentSpeed = norm(velocity_field(:,idx))/norm(original);
    unitDev = velocity_field_dev/norm(velocity_field_dev);
    diffAngle(d) = abs(rad2deg(acos(unitDev'*unitOriginal)));
   

    hold on, quiver(0,0,velocity_field_dev(1), -velocity_field_dev(2),'LineWidth', 3, 'color', on(d,:), 'AutoScaleFactor',1)
    hold on, quiver(0,0,deviation(1), -deviation(2),'LineWidth', 3,'color', off(d,:), 'AutoScaleFactor',1)
end
    hold on
    quiver(0,0, original(1), -original(2),'LineWidth', 3, 'color', 'k', 'AutoScaleFactor',1)
    hold on
    overall = difference(:,1)-difference(:,2);
    quiver(velocity_field_dev(1), -velocity_field_dev(2), overall(1), -overall(2), '--', 'LineWidth', 2,'color', [0 77 128]/255, 'ShowArrowHead', 'off', 'AutoScaleFactor',1)
    
%     quiver(velocity_field_dev(1), -velocity_field_dev(2), difference(1,1), difference(2,1), '--', 'LineWidth', 3,'color', [0 77 128]/255, 'AutoScaleFactor',1)
    title(['Deviations at ' num2str(cond.deviation_depth), ' m'])
    text(.6,.9, ['deviation position = ' mat2str(devPos'), ' deg'], 'Units','normalized')
    text(.6,.85, ['translation = ' mat2str([translate(1), -translate(2), translate(3)]), 'm/s'],'Units','normalized')
    axis equal 
    xlim([0,2.5]), ylim([0,2.5]) 
    yticks([.5 1 1.5 2 2.5])
    
    %for commitee meeting fig
    
    point2 = velocity_field_dev + overall;
    shadex = [0 velocity_field_dev(1) point2(1)];
    shadey = -[0 velocity_field_dev(2) point2(2)];
    hold on, patch(shadex, shadey, 'k', 'LineStyle', 'none', 'FaceAlpha',.2)
    hold on, quiver(velocity_field_dev(1), -velocity_field_dev(2), overall(1), -overall(2), '--', 'LineWidth', 2,'color','k', 'ShowArrowHead', 'off', 'AutoScaleFactor',1)
    
