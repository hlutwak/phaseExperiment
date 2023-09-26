% convert control fig to deg/speed plot

addpath('/Volumes/GoogleDrive/My Drive/opticflow/objectDetection/OpticFlow/phaseExperiment/analysis')
fig = openfig('control_w_conf.fig');

fig = gcf;

axObjs = fig.Children;
dataObjs = axObjs.Children; % figure data

thresholds = findobj(dataObjs, 'Type', 'Line'); %threshold points
xCell = get(thresholds, 'XData');
x = cell2mat(xCell);
yCell = get(thresholds, 'YData');
y = cell2mat(yCell);

baselines = findobj(dataObjs, 'Type', 'Quiver'); %base line velocities tested
uCell = get(baselines, 'UData');
u = cell2mat(uCell);
vCell = get(baselines, 'VData');
v = cell2mat(vCell);

% figure
% for p = 1:length(xCell)
%     hold on
%     scatter(xCell{p}, yCell{p})
% end


basespeeds = vecnorm([u,v], 2, 2); % speeds of baseline velocities tested

threshspeeds = [];
for t = 1:length(xCell) %  speed of thresholds
    threshspeeds(t,:) = vecnorm([xCell{t};yCell{t}]);
end

% figure
% scatter(ones(length(basespeeds),1), basespeeds, 50, 'k', 'filled')
% for t = 1:length(xCell)
%     hold on
%     scatter(ones(1,length(threshspeeds)), threshspeeds(t,:))
% end


basedirections = mod(rad2deg(atan2(v, u)), 360); % directions of baseline velocities
threshdirections = mod(rad2deg(atan2(y, x)), 360); % directions of threshold velocities


%calculate based on difference?
centered_thresh_x = x - u;
centered_thresh_y = y - v;
centered_threshdirections = mod(rad2deg(atan2(centered_thresh_x, centered_thresh_y)), 360); %get directions of thresholds as difference of baseline


figure
scatter(basedirections, 1./basespeeds, 75, 'k', 'filled')
for d = 1:length(xCell)
    hold on
    scatter(threshdirections(d, :), 1./threshspeeds(d,:), 75, 'filled')
end
set(gca,'yscale','log')
% scatter(threshdirections(:), threshspeeds(:), 'k', 'filled');
xlim([0 360])
ylim([0 10])
xlabel('direction')
ylabel('log 1/speed (deg/s)')
