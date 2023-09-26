%% Prediction for detectability across visual field
% gives you a heatmap of how difficult it should be to detect a deviation
% based on the angle between the on and off constraint lines 
% can change translation, relative depth
% addpath(genpath('/Users/hopelutwak/Documents/MATLAB/'))

translate = [ 0, 0.01, 1.];
z0 = .35;
f = .35;
relative_depth = .012;

% visual field
x = linspace(-50, 50, 51);
y = linspace(-30, 30, 31);
[x,y] = meshgrid(x,y);
dots_deg = [x(:)'; y(:)'];

% depths based on relative depth chosen
Z = z0/relative_depth;
dots_m = zeros(3, length(dots_deg));
dots_m(1:2,:) = Z.*tand(dots_deg);
dots_m(3,:) = Z;

% calculate flow field at each point
original_velocity = calculate_cloud_flow(dots_m(3,:), dots_deg, translate, f, z);
dots_m(3,:) = Z/2;
nearer = calculate_cloud_flow(dots_m(3,:), dots_deg, translate, f);
dots_m(3,:) = Z*2;
further = calculate_cloud_flow(dots_m(3,:), dots_deg, translate, f);

diffNear = nearer - original_velocity;
diffFar = further - original_velocity;

% calculate dot product of each original velocity with those of
% nearer/further
dotProdNear = zeros(1,length(original_velocity));
dotProdFar = zeros(1,length(original_velocity));
for ii = 1: length(original_velocity)
    unitOriginal = original_velocity(:,ii)/norm(original_velocity(:,ii));
    unitDiffN = diffNear(:,ii)/norm(diffNear(:,ii));
    unitDiffF = diffFar(:,ii)/norm(diffFar(:,ii));
    
    dotProdNear(ii) = unitOriginal'*unitDiffN;
    dotProdFar(ii) = unitOriginal'*unitDiffF;
end

near = reshape(dotProdNear, 31,51);
far = reshape(dotProdFar, 31, 51);

figure
cmap = cbrewer('div', 'RdBu', 100, 'cubic');
colormap(cmap)

imagesc(near, [-1 1]);  %, [.2 .7]
title(['closer than expected, relative depth = ', num2str(relative_depth), ' original depth = ', num2str(Z)])
xticks([1 26 51])
xticklabels({'-50', '0', '50'})
xlabel('deg')
yticks([1 16 31])
yticklabels({'30', '0', '-30'})
ylabel('deg')
colorbar


% figure
% cmap = cbrewer('div', 'RdBu', 100);
% colormap(cmap)
% imagesc(far);  %, [.2 .7]
% title(['further than expected, relative depth = ', num2str(relative_depth)])
% xticks([1 26 51])
% xticklabels({'-50', '0', '50'})
% xlabel('deg')
% yticks([1 16 31])
% yticklabels({'30', '0', '-30'})
% ylabel('deg')
% colorbar

%%
% dot product between on/off constraint
translate = [ 0, 0.15, .65];
z0 = 3;
f = .35;
middle_dist = 3;
closer_dist = middle_dist*.8;

% visual field
x = linspace(-30, 30, 31);
y = linspace(-24, 24, 25);
[x,y] = meshgrid(x,y);
dots_deg = [x(:)'; y(:)'];

% make filler matrices for calculate_cloud_flow
dots_m = zeros(3, size(dots_deg,2));
dots_m(3,:) = middle_dist;
dots_m_closer = zeros(3,size(dots_deg,2));
dots_m_closer(3,:) = closer_dist;
nDots = length(dots_m);

%calculate velocity fields for each depth
velocity_field_middle = calculate_cloud_flow(dots_m(3,:), dots_deg, translate, f, z0);
velocity_field_closer = calculate_cloud_flow(dots_m_closer(3,:), dots_deg, translate, f, z0);

velocity_field_off = zeros(2,nDots);
%-- get deviation same change in speed and  direction OFF CONSTRAINT-- for
% each location
for ii = 1: nDots
    unitOriginal = velocity_field_middle(:,ii)/norm(velocity_field_middle(:,ii));
    theta = atan(unitOriginal(2)/unitOriginal(1)); %get angle between unit original velocity and x-axis
    R = [cos(theta) sin(theta); -sin(theta) cos(theta)]; % rotates based on tan(theta)
    rotated = R*velocity_field_closer(:,ii); %rotate the velocity along constraint
    rotated(2) = -rotated(2); %reflect over yaxis (rotated direction of original velocity)
    Rback = R'; %rotation matrix to go back
    velocity_field_off(:,ii) = Rback*rotated;
end

%center the contraint segments to compute dot products
centered_constraints_on = velocity_field_closer - velocity_field_middle;
centered_constraints_off = velocity_field_off - velocity_field_middle;
dotProd = zeros(1, nDots);
angle_btwn = zeros(1, nDots);
for ii = 1:nDots
    unit_constraint_on = centered_constraints_on(:,ii)/norm(centered_constraints_on(:,ii));
    unit_constraint_off = centered_constraints_off(:,ii)/norm(centered_constraints_off(:,ii));
    dotProd(ii) = unit_constraint_on'*unit_constraint_off;
    angle_btwn(ii) = rad2deg(acos(dotProd(ii)));
end

% figure
% cmap = cbrewer('div', 'RdBu',100, 'PCHIP');
% colormap(cmap);
% detection = reshape(dotProd, size(x,1), size(x,2));
% imagesc(detection, [-1 1]);
% colorbar;
% set(gca,'XTick', 1:5:x(end)+1, 'XTickLabel', num2str((-x(end): 10: x(end))'));
% set(gca,'YTick', 1:5:y(end)+1, 'YTickLabel', num2str((-y(end): 10: y(end))'));
% 
% title(sprintf('T = [%.2f %.2f %.2f] m/s z0 = %.2fm cloud distance= %dm', translate(1), translate(2), translate(3), z0, middle_dist))
% xlabel('deg')
% ylabel('deg')

figure
cmap = cbrewer('div', 'RdBu',100, 'PCHIP');
colormap(cmap);
detection = reshape(angle_btwn, size(x,1), size(x,2));
imagesc(detection)
colorbar;
set(gca,'XTick', 1:5:x(end)+1, 'XTickLabel', num2str((-x(end): 10: x(end))'));
set(gca,'YTick', 1:5:y(end)+1, 'YTickLabel', num2str((-y(end): 10: y(end))'));
title(sprintf('T = [%.2f %.2f %.2f] m/s z0 = %.2fm cloud distance= %dm', translate(1), translate(2), translate(3), z0, middle_dist))
xlabel('deg')
ylabel('deg')
