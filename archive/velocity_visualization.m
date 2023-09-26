%% Prediction for detectability across visual field
% gives you a heatmap of how difficult it should be to detect a deviation
% based on the dot product of the original velocity, and 
% addpath('/Users/hopelutwak/Documents/opticflow_old/objectDetection/RangeDatabase1080p')
addpath('/Volumes/GoogleDrive/My Drive/opticflow/objectDetection/RangeDatabase1080p')

translate = [0, .5, 1.4];
view_dist                       = .35;                                  % m .57; psychophysics room: .35 .50
screensize                      = [.405 .298];                          % m psychophysics room: [.405 .298], home: [.495 .312]
pixels                          = [1600 1200];                          % psychophysics room: [1600 1200], home: [1920 1200]
frame_rate                      = 60;                                   % Hz 85 , 60 at hhome

view_window                     = round([rad2deg(atan(screensize(1)/2/view_dist)) rad2deg(atan(screensize(2)/2/view_dist))]);  % X,Y centered around fixation [60 46] [54 40];                                                                      % [36 27] for laptop, [55 32] for monitor
scale_factor                    = sqrt((view_window(1)*60*view_window(2)*60)/(pixels(1)*pixels(2)))*1;          % Arcmin/pixel 1.78 2.25 1.92dell: 1680x1050, psychophysics room: 1600x1200                                % Screen frame rate (hz) psychophysics room: 85
scale_factorX = view_window(1)*60/pixels(1);
scale_factorY = view_window(2)*60/pixels(2);
linearize                       = 0;                                    % Use calibrated LUT (do this when available)                    

middle_dist                     = 3;                                   %meters
cloud_dist                      = [middle_dist*.5; middle_dist*1.5];  % m range of distances if depth structure random
bumpiness                       = .5;
devPos                          = [5; 3];
depth_structure                 = 2;                                    % 0 = random cloud, 1 = ground plane, 2 = rangedata
ds                              = [];
    switch depth_structure
        case 0 
        case 1
            
            ds.gaze_angle               = 70;                                   %deg from staring straight at ground
                if (ds.gaze_angle + view_window(2)/2) >=90
                    error('gaze angle too high!')
                end
            translation                 = getT_eyecoords(90-ds.gaze_angle, translation);
            ds.height                   = 1.65;                                 % ds.height in meters
            ds.bumpiness                   = 0.1;                               % random noise in ground plane
            ds.wall                        = 0;                                 % ground plane + wall starting in upper visual field
            ds.planeParallax               = 1;
        case 2
            img_label = 'lRange050.mat';
            load (img_label)
            screen_range_map            = [1.940 1.118];
            screen_range_map_pix        = [1920 1080];
            view_window_map              = rad2deg([atan(screen_range_map(1)/2*1/3), atan(screen_range_map(2)/2*1/3)]); % screen 3 m out
            [frac, axs]                = min(view_window_map./view_window);    % push the scene closer so takes up same amoount of visual angle as experimental screen
            new_dist = screen_range_map(axs)/2*1/tan(deg2rad(view_window(axs)));
            view_window_map = rad2deg([atan(1.94/2*1/new_dist), atan(1.118/2*1/new_dist)]);

            crop_amt = screen_range_map_pix - pixels;              % how much of dispay screen fits on range image *need better solution
            crop_amt(crop_amt<0) = 0;
            crop = rangeMap(1:end - crop_amt(2), 1:end - crop_amt(1), :);
            range_distance = -crop(:,:,1);                                  % swap coordinate system so have positive distances

            range_distance(range_distance == 1) = max(max(range_distance))*1000; % replace all 1s (no laser return - could be sky) with very far distance
            range_distance = range_distance + 3- new_dist;
            
            range_distance_left = range_distance(:,1:pixels(1)/2);
            mirror = [range_distance_left, flip(range_distance_left,2)]; %mirror left side with right so choosing deviation on left and right is symmetric
            range_distance = mirror; 
            pixels = [size(range_distance, 2), size(range_distance, 1)];
            
    end


% visual field
x = linspace(-view_window(1)/2, view_window(1)/2, view_window(1)+1);
y = linspace(-(view_window(2)-1)/2, (view_window(2)-1)/2, view_window(2));
centersurround = [devPos-3, devPos, devPos+3];
% x = centersurround(1,:);
% y = centersurround(2, :);
[xx,yy] = meshgrid(x,y);
dots_deg = [xx(:)'; yy(:)'];
nDots = length(dots_deg);
if depth_structure ==0
    % depths based on relative depth chosen
    Z = middle_dist;
    dots_m = zeros(3, length(dots_deg));
    dots_m(1:2,:) = Z.*tand(dots_deg);
    dots_m(3,:) = Z;
elseif depth_structure == 1
    dots_m = zeros(3, length(dots_deg));
%     dots_m(1:2,:) = view_dist.*tand(dots_deg);
    distances = raytrace2eyeZ(dots_deg, height, gaze_angle);
    
    if wall
        points_upperField = dots_deg(2,:)<0;
        length_groundPlane = height*tan(deg2rad(gaze_angle));
        distances(points_upperField) = raytrace2eyeZ(dots_deg(:,points_upperField), length_groundPlane, -gaze_angle);
    end
    dots_m(3,:) = distances;
    dots_m(1:2,:) = dots_m(3,:).*tand(dots_deg);
    
elseif depth_structure ==2
%     dots_m = zeros(3, length(dots_deg)); 
%     sampley = round(linspace(1,1080, length(y)));
%     samplex = round(linspace(1,1920, length(x)));
%     distances = rangeMap(sampley,samplex,1);  %35 x 20 deg of visual angle
%     depth = abs(distances(:));
%     dots_m(3,:) = depth;
%     dots_m(1:2,:) = dots_m(3,:).*tand(dots_deg);
%     z0 = depth(round(length(depth)/2));
    center  = ceil(pixels/2);
    z0 = range_distance(center(2), center(1));
    devPosPix = [devPos(1)*60*scale_factorX; devPos(2)*60*scale_factorY];
    devPosPix = round(center' + devPosPix);
    middle_dist = range_distance(devPosPix(2), devPosPix(1));
    surround_extent= round([scale_factorX*60*1.1 scale_factorY*60*1.1]); %based on dot density .5 ~2.2 deg diff between center of one dot to next
    surround_window = range_distance(devPosPix(2)-surround_extent(2):min(pixels(2),devPosPix(2)+surround_extent(2)), devPosPix(1)-surround_extent(1):devPosPix(1)+surround_extent(1));

    surround_distance = [min(min(surround_window)), max(max(surround_window))];
    depth_limits = [max(.5, surround_distance(1)), surround_distance(2)]; % take average of depths in window ~2 deg around (based on dot density .5)

    dots_m = zeros(3, length(dots_deg));
    
    dots_deg_pix = [dots_deg(1,:)*60/scale_factorX; dots_deg(2,:)*60/scale_factorY]; %convert to pixels
    dots_deg_pix = round(center'+dots_deg_pix); 

    ok = find(dots_deg_pix(2,:)>0 & dots_deg_pix(2,:)< pixels(2) & dots_deg_pix(1,:)>0 & dots_deg_pix(1,:)<pixels(1)); %get positions which can be found in range map data
    dots_m = dots_m(:,ok);
    dots_deg = dots_deg(:,ok);
    dots_deg_pix = dots_deg_pix(:,ok);

    ind = sub2ind(size(range_distance), dots_deg_pix(2,:), dots_deg_pix(1,:)); %get the indices 
    dots_m(3,:)  = range_distance(ind); %assign depths to distances
    nDots = length(dots_m);

end
% figure
% % dots_m(3,:)+(rand(1, nDots)-.5)*bumpiness;
% plot3(dots_m(1,:), dots_m(2,:), dots_m(3,:), '.')
% hist(dots_m(3,:))

sampling = 1;
% calculate flow field at each point centered at respective point

original_velocity = calculate_cloud_flow(dots_m(3,:), dots_deg, translate, view_dist, z0)+dots_deg;

xMiddleV = [dots_deg(1,1:sampling:end); original_velocity(1,1:sampling:end)];
xMiddleV = xMiddleV(:)';
yMiddleV = [dots_deg(2,1:sampling:end); original_velocity(2,1:sampling:end)];
yMiddleV = yMiddleV(:)';
% figure
% 
% for ii = 1:length(xMiddleV)/2
%     hold on
%     
%     X = xMiddleV(2*ii-1);
%     Y = yMiddleV(2*ii-1);
%     aX = xMiddleV(2*ii);
%     aY = yMiddleV(2*ii);
%     U = xMiddleV(2*ii) - xMiddleV(2*ii-1);
%     V = yMiddleV(2*ii) - yMiddleV(2*ii-1);
% %     X = [xMiddleV(2*ii-1) xMiddleV(2*ii)];
% %     Y = [yMiddleV(2*ii-1) yMiddleV(2*ii)];
% %     plot(X, -Y, 'color',[0.5, .5, .5])
%     quiver(X, -Y, U, -V, 'color',[0, 0, 0], 'linewidth',1, 'AutoScaleFactor', 1)
%       
%     headWidth = 5;
%     headLength = 8;
%     LineLength = 0.08;
%     ah = annotation('arrow',...
%         'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
%     set(ah,'parent',gca);
%     set(ah,'position',[aX -aY LineLength*U -LineLength*V]);
% 
% end
    
if depth_structure == 1 || depth_structure == 2
    cloud_dist = zeros(2, nDots);
    cloud_dist(1,:) = dots_m(3,:)-.5*bumpiness;
    cloud_dist(2,:) = dots_m(3,:)+.5*bumpiness;
end

% calculate flow if all depths were closer than original
dots_m(3,:) = cloud_dist(1,:);
nearer = calculate_cloud_flow(dots_m(3,:), dots_deg, translate, view_dist, z0)+ dots_deg;

% calculate flow if all depths were further than original
dots_m(3,:) = cloud_dist(2,:);
further = calculate_cloud_flow(dots_m(3,:), dots_deg, translate, view_dist, z0)+dots_deg;

diffNear = nearer - original_velocity;
diffFar = further - original_velocity;
%order near velocity and far velocity to connect points

xV = [nearer(1,1:sampling:end); further(1,1:sampling:end)];
xV = xV(:)';
yV = [nearer(2,1:sampling:end); further(2,1:sampling:end)];
yV = yV(:)';

xMiddleV = [dots_deg(1,1:sampling:end); original_velocity(1,1:sampling:end)];
xMiddleV = xMiddleV(:)';
yMiddleV = [dots_deg(2,1:sampling:end); original_velocity(2,1:sampling:end)];
yMiddleV = yMiddleV(:)';

figure
for ii = 1:length(xV)/2
    hold on
    plot([xV(2*ii-1) xV(2*ii)], -[yV(2*ii-1) yV(2*ii)], 'b')
    hold on
    plot([xMiddleV(2*ii-1) xMiddleV(2*ii)], -[yMiddleV(2*ii-1) yMiddleV(2*ii)], 'k', 'Color', [.5 .5 .5])
end


%plot positions where velocity calculated
hold on 
plot(dots_deg(1,1:sampling:end), -dots_deg(2,1:sampling:end), 'k.')
 
%plot deviation positions
hold on
scatter([devPos(1), -devPos(1)], -[devPos(2), devPos(2)], 'r')
axis equal
% xlim([devPos(1)-5, devPos(1)+5])
% ylim([-devPos(2)-5, -devPos(2)+5])
% axis square

% title(sprintf('T = [%.2f %.2f %.2f] m/s z0 = %.2fm cloud distance= %dm', translate(1), translate(2), translate(3), z0, middle_dist))
title(img_label)
xlabel('deg')
ylabel('deg')

%% calculate average of surround

Z = middle_dist;
dots_m = zeros(3, length(dots_deg));
dots_m(1:2,:) = Z.*tand(dots_deg);
dots_m(3,:) = Z;

original_velocity = calculate_cloud_flow(dots_m(3,:), dots_deg, translate, f, z0);
surround_velocity = original_velocity;
surround_velocity(:,4) = [];
surround_mean = mean(surround_velocity,2)


hold on, scatter(surround_mean(1), -surround_mean(2), 50, 'r', 'filled')
%%
figure
colormap gray
imagesc(1./rangeMap(:,:,1))
axis off
%%
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
 
imagesc(far, [-1 1]);  %, [.2 .7]
title(['range of depths ', num2str(cloud_dist),'m'])
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

%% after test run
% look at all velocities on vx vy plot, highlight target velocity red
figure

scatter(velocity_field(1,:), -velocity_field(2,:)), hold on, scatter(velocity_field(1,idx), -velocity_field(2,idx), 'r')
axis equal
if depth_structure ==2
    title(img_label)
else
    title(['deviation position ' num2str(devPos(1)) ',' num2str(devPos(2))])
end

figure
for c = 1:10
    hold on
    scatter(cond(c).steps(1,:), -cond(c).steps(2,:))
    if c == perm(trial)
        scatter(cond(c).steps(1,:), -cond(c).steps(2,:), 'filled')
    end
    axis equal
end
if depth_structure ==2
    title(img_label)
else
    title(['deviation position ' num2str(devPos(1)) ',' num2str(devPos(2))])
end
hold on, c = 1, scatter(cond(c).steps(1,:), -cond(c).steps(2,:),'filled','k')

if depth_structure == 2
    figure, imagesc(range_distance), colorbar, title(img_label)
end

%plot flow map from test run
flowvecs = velocity_field+dots_deg;
figure
for d = 1:length(dots_deg)
    hold on
    plot([dots_deg(1,d) flowvecs(1,d)], -[dots_deg(2,d) flowvecs(2,d)], 'k', 'Color', [.5 .5 .5])
end
hold on
plot(dots_deg(1,:), -dots_deg(2,:), 'k.')
axis equal
 
%plot deviation positions
hold on
scatter([devPos(1), -devPos(1)], -[devPos(2), devPos(2)], 'r')


