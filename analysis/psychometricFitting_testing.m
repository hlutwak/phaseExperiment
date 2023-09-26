%% Psychometric fitting
addpath('/Users/hopelutwak/Documents/MATLAB/psignifit')
%plot staircases
c = jet(length(cond));
c = c*.95;
angles = [0 45 90 135 180 225 270 315];
num_staircases = 2;
num_t = 1;
steps = 1:length(cond(1).steps);

figure

for ii = 1:length(cond)/num_staircases
     subplot(num_t,length(angles),ii)
     plot(-cond(2*ii-1).step_history, '-o', 'color', colors(2*ii-1,:), 'LineWidth',2, 'DisplayName', num2str(cond(2*ii-1).rotation))
     hold on
     plot(-cond(2*ii).step_history, '-o', 'color', colors(2*ii,:), 'LineWidth',2, 'DisplayName', num2str(cond(2*ii).rotation))
     title(num2str(cond(2*ii).rotation))
     ylim([-10 0])
end

% plot in terms of magnitude delta
% deltanorm = (vecnorm(cond(1).steps - cond(1).steps(:,end)));
% figure
% for ii = 1:length(cond)/num_staircases
%      subplot(num_t,length(angles),ii)
%      plot(-cond(2*ii-1).step_history, '-o', 'color', c(2*ii-1,:), 'LineWidth',2, 'DisplayName', num2str(cond(2*ii-1).rotation))
%      hold on
%      plot(-cond(2*ii).step_history, '-o', 'color', c(2*ii,:), 'LineWidth',2, 'DisplayName', num2str(cond(2*ii).rotation))
%      title(num2str(cond(2*ii).rotation))
% end


% plot resp in terms of weber fractions
% figure
% for ii = 1:length(angles)
%      subplot(2,4,ii)
%      plot((abs(vecnorm(cond(2*ii-1).dev_history) - norm(cond(2*ii-1).steps(:,end))))/ norm(cond(2*ii-1).steps(:,end)), '-o', 'color', c(2*ii-1,:), 'LineWidth',2, 'DisplayName', num2str(cond(2*ii-1).rotation))
%      hold on
%      plot((abs(vecnorm(cond(2*ii).dev_history) - norm(cond(2*ii).steps(:,end))))/ norm(cond(2*ii).steps(:,end)), '-o', 'color', c(2*ii,:), 'LineWidth',2, 'DisplayName', num2str(cond(2*ii).rotation))
%      title(num2str(cond(2*ii).rotation))
%      ylim([0 1])
% end
% 
% % calculate weber fraction
% weberFrac = zeros(length(angles), length(cond(1).steps));
% for ii = 1:length(angles)
%     weberFrac(ii,:) = (abs(vecnorm(cond(2*ii).steps) - norm(cond(2*ii).steps(:,end))))/ norm(cond(2*ii).steps(:,end));
% end
% % 
% % get velocity diff for each condition
% delta = cell(length(angles), 1);
% for ii = 1:length(angles)
%     delta{ii} = cond(2*ii).steps - cond(2*ii).steps(:,end);
% end


figure
% plot velocitities tested
for ii = 1:length(cond)/num_staircases
    hold on, scatter(-cond(ii*num_staircases).steps(1,:), -cond(ii*num_staircases).steps(2,:), 50, [.5 .5 .5], 'filled')
end
hold on, quiver(0, 0, -cond(1).steps(1,end), -cond(1).steps(2,end), 'Color', 'k','LineWidth', 2,'AutoScaleFactor',1)
axis equal


%% percent correct
% load data
% myFolder = '/Volumes/GoogleDrive/My Drive/opticflow/objectDetection/OpticFlow/phaseExperiment/data/'

% because I messed up real bad, recalculate num correct 
%*fixed bug in experiment code 11/24
%for now get pcorrect for everything but last trial
% for ii = 1:length(cond)
%     cond(ii).resp_history = [];
%     cond(ii).resp_history = diff(cond(ii).step_history)> -1;
%     cond(ii).step_history(:,end) = [];
% end

deltanorm = (vecnorm(cond(1).steps - cond(1).steps(:,end)));
% get percent correct for each step
nCorrect = zeros(length(angles)*num_t,length(cond(1).steps));
nTrials = zeros(length(angles)*num_t,length(cond(1).steps));


% for ii = 1:length(cond)/num_staircases
%     for step = 1:length(cond(1).steps)
%         
%            idx = find(cond(2*ii-1).step_history ==step);
%            numTrials = length(idx);
%            ncorrect = sum(cond(2*ii-1).resp_history(idx));
%            
%            idx = find(cond(2*ii).step_history ==step);
%            numTrials2 = length(idx);
%            ncorrect2 = sum(cond(2*ii).resp_history(idx));
%            
%            nCorrect(ii,step) = nCorrect(ii,step)+ncorrect+ncorrect2;
%            nTrials(ii,step) = nTrials(ii,step)+numTrials + numTrials2;
%     end
% end

%for one staircase in a block
for ii = 1:length(cond)
    for step = 1:length(cond(1).steps)
        
           idx = find(cond(ii).step_history ==step);
           numTrials = length(idx);
           ncorrect = sum(cond(2*ii-1).resp_history(idx));
           
           nCorrect(ii,step) = nCorrect(ii,step)+ncorrect;
           nTrials(ii,step) = nTrials(ii,step)+numTrials;
    end
end

nTrials(nTrials==0) = NaN;
pCorrect = nCorrect./nTrials;


for ii = 1:length(cond)/num_staircases
     subplot(num_t,length(angles),ii)
     plot(-cond(2*ii-1).step_history, '-o', 'color', colors(2*ii-1,:), 'LineWidth',2, 'DisplayName', num2str(cond(2*ii-1).rotation))
     hold on
     plot(-cond(2*ii).step_history, '-o', 'color', colors(2*ii,:), 'LineWidth',2, 'DisplayName', num2str(cond(2*ii).rotation))
     title(num2str(cond(2*ii).rotation))
     ylim([-10 0])
end

figure
for ii = 1:length(cond)/num_staircases
    subplot(num_t, length(angles),ii)
    plot(deltanorm,100*pCorrect(ii,:),'-','MarkerFaceColor','b');
%     plot(-steps,100*pCorrect(ii,:),'-','MarkerFaceColor','b');

     %loop through each intensity so each data point can have it's own size.
    for jj= steps
        if ~isnan(pCorrect(ii,jj))
            sz = nTrials(ii,jj)+1;
            hold on
            plot(deltanorm(jj),100*pCorrect(ii,jj),'ko-','MarkerFaceColor','b','MarkerSize',sz);
%             plot(-jj,100*pCorrect(ii,jj),'ko-','MarkerFaceColor','b','MarkerSize',sz);
        end
    end
    set(gca,'YLim',[0,100]);
end

set(gca,'XTick',steps);
set(gca,'YLim',[0,100]);
% xlabel('staircase');
ylabel('Percent Correct');


%% fit psychometric curve
angle = 315; %which angle
angles = [0 45 90 135 180 225 270 315];
idx = find(angles == angle);
% put responses and intensities in terms of delta into vectors
responses = [cond(num_staircases*idx-1).resp_history cond(num_staircases*idx).resp_history];
intensities = [deltanorm(cond(num_staircases*idx-1).step_history) deltanorm(cond(num_staircases*idx-1).step_history)];

figure
for jj= steps
    if ~isnan(pCorrect(idx,jj))
        sz = nTrials(idx,jj)+1;
        hold on
        plot(deltanorm(jj),100*pCorrect(idx,jj),'ko-','MarkerFaceColor','b','MarkerSize',sz);
    end
end
set(gca,'YLim',[0,100]);

t = .7;
s = 4;
deltas = linspace(min(deltanorm), max(deltanorm), 100);
% get log likelihood for our guess
negloglikelihood = negloglike(t, s, intensities, responses);

y = Weibull(t, s, deltas);
hold on
plot(deltas, y*100, 'r-', 'LineWidth', 2, 'DisplayName',['guess, t = ' num2str(t), ' negloglike = ' num2str(negloglikelihood)])

% fminsearch
f = @(x) negloglike(x(1),x(2),intensities,responses);
x0(1) = 1;
x0(2) = 1;
[x_opt, fval] = fminsearch(f, x0);

y = Weibull(x_opt(1), x_opt(2), deltas);
hold on
plot(deltas, y*100, 'g-', 'LineWidth', 2, 'DisplayName',['search, t = ' num2str(x_opt(1)), ' negloglike = ' num2str(fval)])
legend
% contour
figure
f = @(x1, x2) negloglike(x1,x2,intensities,responses);
fcontour(f, [0 1.5 0 6]), colorbar


%% fminsearch over conditions

% setup
c = jet(length(cond));
c = c*.95;
angles = [0 45 90 135 180 225 270 315];
num_staircases = 2;
num_t = 1;
steps = 1:length(cond(1).steps);
% convert steps into measure of intensity (magnitude of difference vector)
deltanorm = (vecnorm(cond(1).steps - cond(1).steps(:,end)));

% load data

% IF BEFORE 11/24
% because I messed up, recalculate num correct 
% get pcorrect for everything but last trial
for ii = 1:length(cond)
    cond(ii).resp_history = [];
    cond(ii).resp_history = diff(cond(ii).step_history)> -1;
    cond(ii).step_history(:,end) = [];
end

% put responses and intensities in terms of delta into cell array results
responses = [cond(1).resp_history cond(2).resp_history];
intensities = [deltanorm(cond(1).step_history) deltanorm(cond(2).step_history)];

f = @(x) negloglike(x(1),x(2),intensities,responses);

x0(1) = 1;
x0(2) = 2;
x_opt = fminsearch(f, x0)


%% psignifit

addpath('/Users/hopelutwak/Documents/MATLAB/psignifit')
%/Volumes/GoogleDrive/My Drive/opticflow/objectDetection/OpticFlow/phaseExperiment

% IF BEFORE 11/24
% because I messed up, recalculate num correct 
% get pcorrect for everything but last trial
% for ii = 1:length(cond)
%     cond(ii).resp_history = [];
%     cond(ii).resp_history = diff(cond(ii).step_history)> -1;
%     cond(ii).step_history(:,end) = [];
% end

c = jet(length(cond));
c = c*.95;
angles = [0 45 90 135 180 225 270 315 0 180];
num_staircases = 2;
num_t = 1;
steps = 1:length(cond(1).steps);
% convert steps into measure of intensity (magnitude of difference vector)
deltanorm = (vecnorm(cond(1).steps - cond(1).steps(:,end)));


% plot staircases
figure
for ii = 1:length(cond)/num_staircases
     subplot(num_t,length(angles),ii)
     plot(-cond(2*ii-1).step_history, '-o', 'color', c(2*ii-1,:), 'LineWidth',2, 'DisplayName', num2str(cond(2*ii-1).rotation))
     hold on
     plot(-cond(2*ii).step_history, '-o', 'color', c(2*ii,:), 'LineWidth',2, 'DisplayName', num2str(cond(2*ii).rotation))
     title(num2str(cond(2*ii).rotation))
     ylim([-10 0])
end

% plot velocitities tested
figure
for ii = 1:length(cond)/num_staircases
    hold on, scatter(cond(ii*num_staircases).steps(1,:), -cond(ii*num_staircases).steps(2,:), 50, c(2*ii-1,:), 'filled')
end
hold on, quiver(0, 0, cond(1).steps(1,end), -cond(1).steps(2,end), 'Color', 'k','LineWidth', 2,'AutoScaleFactor',1)
axis equal


%
% get n correct for each stim
nCorrect = zeros(length(angles)*num_t,length(cond(1).steps));
nTrials = zeros(length(angles)*num_t,length(cond(1).steps));

if num_staircases == 2
    for ii = 1:length(cond)/num_staircases
        for step = 1:length(cond(1).steps)
            
            idx = find(cond(2*ii-1).step_history ==step);
            numTrials = length(idx);
            ncorrect = sum(cond(2*ii-1).resp_history(idx));
            
            idx = find(cond(2*ii).step_history ==step);
            numTrials2 = length(idx);
            ncorrect2 = sum(cond(2*ii).resp_history(idx));
            
            nCorrect(ii,step) = nCorrect(ii,step)+ncorrect+ncorrect2;
            nTrials(ii,step) = nTrials(ii,step)+numTrials + numTrials2;
        end
    end
    
else
    for ii = 1:length(cond)
        for step = 1:length(cond(1).steps)
            
            idx = find(cond(ii).step_history ==step);
            numTrials = length(idx);
            ncorrect = sum(cond(2*ii-1).resp_history(idx));
            
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
%loop through stim conditions and get threshold and plot curves
figure
thresholds = zeros(1,length(angles)*num_t);
CI = zeros(2, length(angles)*num_t);
for ii = 1:length(angles)*num_t
    data = [flip(steps)' nCorrect(ii,:)' nTrials(ii,:)'];
    data(data(:,end)==0, :) = [];
    result = psignifit(data,options);
    thresholds(ii)= exp(result.Fit(1));
    CI(:,ii) =exp(result.conf_Intervals(1,:,1))'; %get 95% CI (first layer) for threshold (first row)
    
    subplot(num_t,length(angles),ii)
    
    plotPsych(result);
    ylim([0 1])
end


% plot thresholds in velocity space

% convert threshold and CI to velocity
velocity_thresh = zeros(2,length(angles)*num_t);
velocity_CI_low = zeros(2, length(angles)*num_t);
velocity_CI_high = zeros(2, length(angles)*num_t);

for ii = 1:length(angles)*num_t
    segment = cond(num_staircases*ii).steps(:,1) - cond(num_staircases*ii).steps(:,end);
    velocity_thresh(:,ii) = cond(num_staircases*ii).steps(:,end)+0.1*thresholds(ii)*segment;
    velocity_CI_low(:,ii) = cond(num_staircases*ii).steps(:,end)+0.1*CI(1,ii)*segment;
    velocity_CI_high(:,ii) = cond(num_staircases*ii).steps(:,end)+0.1*CI(2,ii)*segment;
end
figure
hold on, fill([velocity_CI_high(1,1:length(angles)) velocity_CI_high(1,1)],-[velocity_CI_high(2,1:length(angles)) velocity_CI_high(2,1)], [.5, .5, 1], 'LineStyle','none', 'FaceAlpha',0.5)
hold on, fill([velocity_CI_low(1,1:length(angles)) velocity_CI_low(1,1)],-[velocity_CI_low(2,1:length(angles)) velocity_CI_low(2,1)], 'w', 'LineStyle','none')

hold on, plot([velocity_thresh(1,1:length(angles)) velocity_thresh(1,1)],-[velocity_thresh(2,1:length(angles)) velocity_thresh(2,1)], 'color',[.2, .2, .7], 'linewidth', 2)
hold on, quiver(0, 0, cond(1).steps(1,end), -cond(1).steps(2,end), 'Color', 'k','LineWidth', 2,'AutoScaleFactor',1)
axis equal

if num_t ==2
    figure
    hold on, plot([velocity_thresh(1,length(angles)+1:end) velocity_thresh(1,length(angles)+1)],-[velocity_thresh(2,length(angles)+1:end) velocity_thresh(2,length(angles)+1)], 'color',[.7, .5, .5], 'linewidth', 2)
    hold on, quiver(0, 0, cond(1).steps(1,end), -cond(1).steps(2,end), 'Color', 'k','LineWidth', 2,'AutoScaleFactor',1)
    axis equal
end

%%
figure, 
hold on, scatter(velocity_thresh(1,:), -velocity_thresh(2,:), 50, [.7, .7, .7], 'filled')
hold on, plot(cond(1).velocity_range(1,:), -cond(1).velocity_range(2,:), 'color', [0.5 0.5 0.5])
axis equal

%% Dec 2020

addpath('/Users/hopelutwak/Documents/MATLAB/psignifit')
%/Volumes/GoogleDrive/My Drive/opticflow/objectDetection/OpticFlow/phaseExperiment
% 
% load('data/2020-Dec-05_HL_tup')
left = 0;

num_staircases = 2;
num_t = 1;
steps = 1:length(cond(1).steps);
% convert steps into measure of intensity (magnitude of difference vector)
deltanorm = (vecnorm(cond(1).steps - cond(1).steps(:,end)));

colors = jet(length(cond));
colors = colors*.95;

% plot velocitities tested
figure(1)
for ii = 1:length(cond)/num_staircases
    hold on, scatter(cond(ii*num_staircases).steps(1,:), -cond(ii*num_staircases).steps(2,:), 50, colors(num_staircases*ii,:), 'filled')
end
hold on, quiver(0, 0, cond(1).steps(1,end), -cond(1).steps(2,end), 'Color', 'k','LineWidth', 2,'AutoScaleFactor',1)
axis equal

% get angles, and reorder conditions so radially ordered
angles = zeros(1,length(cond));
for c = 1:length(cond)
    % get dir dev in same dir as base velocity
    base = cond(1).steps(:,1) - cond(1).steps(:,end);
    unit_base = [base/norm(base);0];
    deviation = cond(c).steps(:,1) - cond(c).steps(:,end);
    unit_deviation = [deviation/norm(deviation);0];
    cosTheta = unit_base'*unit_deviation;
    crossP = cross(unit_base, unit_deviation);
        if left
            crossP = -crossP;
        end
    sinTheta = crossP(end);

    angles(c) = mod(-sign(sinTheta)*rad2deg(acos(cosTheta)), 360);
end
[sorted, angleidx] = sort(angles);
n_angles = length(unique(angles));
sorted_unique = unique(sorted);


% plot staircases
figure(2)
for ii = 1:length(cond)/num_staircases
     subplot(num_t,length(angles)/num_staircases,ii)
     anchor_idx = ismember(1:length(cond(1).step_history), cond(num_staircases*ii).anchors);
     hold on
     plot(-cond(num_staircases*ii).step_history(~anchor_idx), '-o', 'color', colors(num_staircases*ii,:), 'LineWidth',2, 'DisplayName', num2str(cond(num_staircases*ii).rotation))
     title(num2str(angles(num_staircases*ii)))
     if num_staircases >1
        anchor_idx = ismember(1:length(cond(1).step_history), cond(num_staircases*ii-1).anchors);
        hold on, plot(-cond(num_staircases*ii-1).step_history(~anchor_idx), '-o', 'color', colors(num_staircases*ii-1,:), 'LineWidth',2, 'DisplayName', num2str(cond(num_staircases*ii-1).rotation))
     
     end
     ylim([-10 0])
end
%
% get n correct for each stim
nCorrect = zeros(length(unique(angles))*num_t,length(cond(1).steps));
nTrials = zeros(length(unique(angles))*num_t,length(cond(1).steps));

sample_angle = angleidx(num_staircases:num_staircases:end);

%save in correct angle order
if num_staircases == 2
    for ii = 1:length(cond)/num_staircases
        for step = 1:length(cond(1).steps)
            
            idx = find(cond(sample_angle(ii)-1).step_history ==step);
            numTrials = length(idx);
            ncorrect = sum(cond(sample_angle(ii)-1).resp_history(idx));
            
            idx = find(cond(sample_angle(ii)).step_history ==step);
            numTrials2 = length(idx);
            ncorrect2 = sum(cond(sample_angle(ii)).resp_history(idx));
            
            nCorrect(ii,step) = nCorrect(ii,step)+ncorrect+ncorrect2;
            nTrials(ii,step) = nTrials(ii,step)+numTrials + numTrials2;
        end
    end
    
else
    for ii = 1:length(cond)
        for step = 1:length(cond(1).steps)
            
            idx = find(cond(sample_angle(ii)).step_history ==step);
            numTrials = length(idx);
            ncorrect = sum(cond(sample_angle(ii)).resp_history(idx));
            
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

%loop through stim conditions and get threshold and plot curves
figure(3)
thresholds = zeros(1,length(angles)/num_staircases*num_t);
CI = zeros(2, length(angles)/num_staircases*num_t);
for ii = 1:length(unique(angles))*num_t
    data = [flip(steps)' nCorrect(ii,:)' nTrials(ii,:)'];
    data(data(:,end)==0, :) = [];
    result = psignifit(data,options);
    thresholds(ii)= exp(result.Fit(1));
    CI(:,ii) =exp(result.conf_Intervals(1,:,1))'; %get CI for threshold (first row) at 95% (first layer)
    
    subplot(num_t,length(unique(angles)),ii)
    
    plotPsych(result);
    ylim([0 1])
    title(num2str(sorted_unique(ii)))
end


% plot thresholds in velocity space

% convert threshold and CI to velocity
velocity_thresh = zeros(2,n_angles*num_t);
velocity_CI_low = zeros(2, n_angles*num_t);
velocity_CI_high = zeros(2, n_angles*num_t);


%save in angle order
for ii = 1:length(cond)/num_staircases
    segment = cond(sample_angle(ii)).steps(:,1) - cond(sample_angle(ii)).steps(:,end);
    velocity_thresh(:,ii) = cond(sample_angle(ii)).steps(:,end)+0.1*thresholds(ii)*segment;
    velocity_CI_low(:,ii) = cond(sample_angle(ii)).steps(:,end)+0.1*CI(1,ii)*segment;
    velocity_CI_high(:,ii) = cond(sample_angle(ii)).steps(:,end)+0.1*CI(2,ii)*segment;
end

%plot boudnary with confidence intervals
colors = [[.5 .5 .5];[.8 .8 .8];[.2 .2 .7];[.7 .7 1]];
figure

% plot CIs
hold on, fill([velocity_CI_high(1,:) velocity_CI_high(1,1)],-[velocity_CI_high(2,:) velocity_CI_high(2,1)], colors(4,:), 'LineStyle','none', 'FaceAlpha',0.5)
hold on, fill([velocity_CI_low(1,:) velocity_CI_low(1,1)],-[velocity_CI_low(2,:) velocity_CI_low(2,1)], 'w', 'LineStyle','none')
% plot thresholds
hold on, plot([velocity_thresh(1,:) velocity_thresh(1,1)],-[velocity_thresh(2,:) velocity_thresh(2,1)], 'color',colors(3,:), 'linewidth', 2)
% plot base velocity
hold on, quiver(0, 0, cond(1).steps(1,end), -cond(1).steps(2,end), 'Color', 'k','LineWidth', 2,'AutoScaleFactor',1)
% plot constraint ine
hold on, plot(cond(1).velocity_range(1,:), -cond(1).velocity_range(2,:), 'color', colors(1,:), 'linewidth', 2)
axis equal

if num_t ==2
    figure
    hold on, plot([velocity_thresh(1,n_angles+1:end) velocity_thresh(1,n_angles+1)],-[velocity_thresh(2,n_angles+1:end) velocity_thresh(2,n_angles+1)], 'color',[.7, .5, .5], 'linewidth', 2)
    hold on, quiver(0, 0, cond(1).steps(1,end), -cond(1).steps(2,end), 'Color', 'k','LineWidth', 2,'AutoScaleFactor',1)
    axis equal
end

