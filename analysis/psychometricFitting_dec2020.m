% Dec 2020
% add psignifit toolbox
addpath('/Users/hopelutwak/Documents/MATLAB/psignifit')

% assign data folder
dataFolder = '/Volumes/GoogleDrive/My Drive/opticflow/objectDetection/OpticFlow/phaseExperiment/data';

% names of files
S = dir(fullfile(dataFolder,'*.mat'));

% which subject data to analyze
subject = 'HL';
stim = 'C_natural'; % 'T1' 'T2'

%load appropriate files
count = 1;
for f = 1:length(S)
    subj = contains(S(f).name,subject);
    sti = contains(S(f).name,stim);
    if subj && sti
        file{count} = load(S(f).name);
        count = count+1;
    end
end

if contains(stim, 'C')
    control =1;
else 
    control = 0;
end
% scramble = 0;
t = 1;
num_t = 1;
steps = 1:size(cond(1).steps,2);
n_blocks = 4;
n_conditions = length(cond)/n_blocks;
n_steps = length(steps);
% convert steps into measure of intensity (magnitude of difference vector)
deltanorm = (vecnorm(cond(1).steps(:,:,1) - cond(1).steps(:,end,1)));

colors = jet(length(cond)/n_blocks);
colors = colors*.95;

% plot velocitities tested
figure(4)
for ii = 1:length(cond)
    hold on, scatter(cond(ii).steps(1,:), -cond(ii).steps(2,:), 50, colors(ii - (cond(ii).block-1)*n_steps,:), 'filled')
end
hold on, quiver(0, 0, cond(1).steps(1,end,1), -cond(1).steps(2,end,1), 'Color', 'k','LineWidth', 2,'AutoScaleFactor',1)
axis equal

% get angles, and reorder conditions so radially ordered
deltaAngles = zeros(1,length(cond));

for c = 1:length(cond)
    % get dir dev in same dir as base velocity
    base = cond(11).steps(:,1,1) - cond(11).steps(:,end,1);
    unit_base = [base/norm(base);0];
    deviation = cond(c).steps(:,1,1) - cond(c).steps(:,end,1);
    unit_deviation = [deviation/norm(deviation);0];

        if control && mod(cond(c).block, 2) && ~scramble
            unit_deviation(1) = -unit_deviation(1);
        end
    cosTheta = unit_base'*unit_deviation;

    crossP = cross(unit_base, unit_deviation);
    sinTheta = crossP(end);
    deltaAngles(c) = mod(-sign(sinTheta)*rad2deg(real(acos(cosTheta))), 360);

end
[sorted_cond, angleidx] = sort(deltaAngles);
n_angles = length(unique(deltaAngles));
sorted_unique_deltaAngles = unique(sorted_cond);



% plot staircases
figure

for ii = 1:n_conditions %n_conditions = n_angles
     % which staircases are have this angle
     cond_idx = angleidx(ii*n_blocks- (n_blocks-1):ii*n_blocks);
     
     subplot(num_t,n_conditions,ii)
     for jj = cond_idx % go through these staircases and plot
         anchor_idx = ismember(1:length(cond(1).step_history), cond(jj).anchors);
         hold on
         plot(-cond(jj).step_history(~anchor_idx), '-o', 'color', colors(ii,:), 'LineWidth',2)
     end
     title(num2str([sorted_unique_deltaAngles(ii)]));
     ylim([-10 0])
end

%
% get n correct for each stim
nCorrect = zeros(n_conditions,size(cond(1).steps,2));
nTrials = zeros(n_conditions,size(cond(1).steps,2));

%save in correct angle order

for ii = 1:n_conditions
    % which staircases are have this angle
    cond_idx = angleidx(ii*n_blocks- (n_blocks-1):ii*n_blocks);
    
    for jj = cond_idx % go through this staircase and count num correct
        for step = 1:size(cond(1).steps,2)
            idx = find(cond(jj).step_history ==step);
            numTrials = length(idx);
            ncorrect = sum(cond(jj).resp_history(idx));

            nCorrect(ii,step) = nCorrect(ii,step)+ncorrect;
            nTrials(ii,step) = nTrials(ii,step)+numTrials;
        end
    end
end


% set options for psychometric functions
options             = struct;   % initialize as an empty struct
options.sigmoidName = 'weibull';   
options.expType     = '2AFC';   % choose 2-AFC as the paradigm of the experiment
                                % this sets the guessing rate to .5 and
                                % fits the rest of the parameters
options.fixedPars = NaN(5,1);                                
options.fixedPars(5) = 0;       % fix eta (dispersion) to zero

%loop through stim conditions and get threshold and plot curves
figure
thresholds = zeros(1,n_conditions);
CI = zeros(2, n_conditions);
for ii = 1:n_conditions
    data = [flip(steps)' nCorrect(ii,:)' nTrials(ii,:)'];
    data(data(:,end)==0, :) = [];
    result = psignifit(data,options);
    thresholds(ii)= exp(result.Fit(1));
    CI(:,ii) =exp(result.conf_Intervals(1,:,1))'; %get CI for threshold (first row) at 95% (first layer)
    
    subplot(num_t,n_conditions,ii)
    
    plotPsych(result);
    ylim([0 1])
    title(num2str(sorted_unique_deltaAngles(ii)))
end


% plot thresholds in velocity space

% convert threshold and CI to velocity
velocity_thresh = zeros(2,n_angles*num_t);
velocity_CI_low = zeros(2, n_angles*num_t);
velocity_CI_high = zeros(2, n_angles*num_t);

%save in angle order
for ii = 1:n_conditions
    % which staircases are have this angle
    cond_idx = angleidx(ii*n_blocks- (n_blocks-1):ii*n_blocks);
    
    segment = cond(cond_idx(2)).steps(:,1,1) - cond(cond_idx(2)).steps(:,end,1);
    
    velocity_thresh(:,ii) = cond(cond_idx(2)).steps(:,end,1)+0.1*thresholds(ii)*segment;
    velocity_CI_low(:,ii) = cond(cond_idx(2)).steps(:,end,1)+0.1*CI(1,ii)*segment;
    velocity_CI_high(:,ii) = cond(cond_idx(2)).steps(:,end,1)+0.1*CI(2,ii)*segment;
    
end

%plot boudnary with confidence intervals
figure

hold on
if control
    colors = [[.5 .5 .5];[.8 .8 .8]];
elseif t ==2
    colors = [[.2 .7 .2];[.8 1 .8]];
else
    colors = [[.2 .2 .7];[.7 .7 1]];
end

% plot CIs
hold on, fill([velocity_CI_high(1,:) velocity_CI_high(1,1)],-[velocity_CI_high(2,:) velocity_CI_high(2,1)], colors(2,:), 'LineStyle','none', 'FaceAlpha',0.5)
hold on, fill([velocity_CI_low(1,:) velocity_CI_low(1,1)],-[velocity_CI_low(2,:) velocity_CI_low(2,1)], 'w', 'LineStyle','none')
% plot thresholds
hold on, plot([velocity_thresh(1,:) velocity_thresh(1,1)],-[velocity_thresh(2,:) velocity_thresh(2,1)], 'color',colors(1,:), 'linewidth', 2)
% plot base velocity
hold on, quiver(0, 0, cond(cond_idx(2)).steps(1,end,1), -cond(cond_idx(2)).steps(2,end,1), 'Color', 'k','LineWidth', 2,'AutoScaleFactor',1)
% % plot constraint line
if ~control
    hold on, plot(cond(cond_idx(2)).velocity_range(1,:,1), -cond(cond_idx(2)).velocity_range(2,:,1), 'color', [.5 .5 .5], 'linewidth', 2)
end

axis equal

if num_t ==2
    figure
    hold on, plot([velocity_thresh(1,n_angles+1:end) velocity_thresh(1,n_angles+1)],-[velocity_thresh(2,n_angles+1:end) velocity_thresh(2,n_angles+1)], 'color',[.7, .5, .5], 'linewidth', 2)
    hold on, quiver(0, 0, cond(cond_idx(2)).steps(1,end,1), -cond(cond_idx(2)).steps(2,end,1), 'Color', 'k','LineWidth', 2,'AutoScaleFactor',1)
    axis equal
end
% xlim([-.2, .8])
% ylim([-.4,.2])
% 
% % plot magnitude distance from base velocity
% % get angle of constraint line deviations
% for c = 1:8
%     theta(c) = cond(c).rotation;
% end
% devAngles = setdiff(round(sorted_unique_deltaAngles,5), theta);

% distance = vecnorm(velocity_thresh - cond(cond_idx(2)).steps(1,end,1));
% figure
% hax = axes;
% hold on
% plot(sorted_unique_deltaAngles, distance,  'color',colors(1,:),'LineWidth', 2)
% hold on
% xline(devAngles(1), 'r:', 'LineWidth', 2)
% hold on
% xline(devAngles(2), 'r:', 'LineWidth', 2)
% hold on
% xline(180, 'k--', 'LineWidth', 2)
% xlabel('angle between axis of tested velocities to base velocity')
% ylabel('magnitude difference between threshold velocity and base velocity')

% % plot fractional speed vs angle diff vs threshold magnitude
% angles = zeros(length(velocity_thresh), 1);
% speeds = zeros(length(velocity_thresh),1);
% for t = 1:length(velocity_thresh)
%     
%     unit_base_velocity = cond(cond_idx(2)).steps(:,end,1)/norm(cond(cond_idx(2)).steps(:,end,1));
%     unit_velocity = velocity_thresh(:,t)/norm(velocity_thresh(:,t));
%     angles(t) = rad2deg(acos(unit_base_velocity'*unit_velocity));
%     
%     speeds(t) = norm(velocity_thresh(:,t))/norm(cond(cond_idx(2)).steps(:,end,1));
% end

% figure
% hold on
% scatter3(angles, speeds, distance,50,colors(1,:), 'filled')
% figure
% hold on
% plot3(angles, speeds, distance, 'LineWidth', 2, 'color', colors(1,:))
% xlabel('angle')
% ylabel('speed')
% zlabel('magnitude threshold')

% figure
% plot(angles, distance, 'color', colors(1,:), 'LineWidth', 2)
% xlabel('angle')
% ylabel('thresh')
% 
% figure
% plot(speeds, distance, 'color', colors(1,:), 'LineWidth', 2)
% xlabel('speed')
% ylabel('thresh')

% plot change in speed (fraction of original) vs change in angle

    