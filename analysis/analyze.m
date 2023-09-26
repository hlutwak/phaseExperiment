%% Analyse data
%%
load ('./data/hl_data.mat')
% import dms file
data = table2array(hl_data);

 
upOnidx = (data(:,3)< 0 & data(:,5) == 0);
upOffidx = (data(:,3)< 0 & data(:,5) == 1);
downOnidx = (data(:,3)> 0 & data(:,5) == 0);
downOffidx = (data(:,3)> 0 & data(:,5) == 1);
 
upOn = data(upOnidx,:);
upOff = data(upOffidx,:);
downOn = data(downOnidx,:);
downOff = data(downOffidx,:);

figure
plot(upOn(:,7))
hold on
plot(upOff(:,7))
hold on
plot(downOn(:,7))
hold on
plot(downOff(:,7))
legend('off down', 'on down', 'off up','on up')

figure
plot(upOn(:,8))
hold on
plot(upOff(:,8))
hold on
plot(downOn(:,8))
hold on
plot(downOff(:,8))
legend('off down', 'on down', 'off up','on up')

%%

c = flipud(jet(length(cond)));
c = c*.9;
conditions = {'down off','down off', 'down on', 'down on','up off', 'up off', 'up on','up on'};
% conditions = {'right off', 'right off', 'right on','right on', 'left off','left off', 'left on', 'left on'};
% cidx = [1 2 5 6 3 4 7 8];
cidx = [1 2 3 4 5 6 7 8];
%depth
Legend = cell(length(cond), 1);
figure 
    for ii = 1:length(cond)
        hold on
        if ii == 3 || ii == 4 || ii == 7|| ii == 8
            plot(cond(ii).dev_depth_history, '-.o', 'LineWidth',2,  'color', c(cidx(ii),:))
        else
            plot(cond(ii).dev_depth_history, '-o', 'LineWidth',2,  'color', c(cidx(ii),:))
        end
        Legend{ii} = [conditions{ii} ' ' num2str(mean(cond(ii).dev_depth_history(15:end)))];
    end
legend('Location', 'eastoutside')
title('depth')
 
 
%angle
figure 
    for ii = 1:length(cond)
        hold on
        if ii == 3 || ii == 4 || ii == 7|| ii == 8
            plot(cond(ii).dev_angle, '-.o', 'LineWidth',2,  'color', c(cidx(ii),:))
        else
            plot(cond(ii).dev_angle, '-o', 'LineWidth',2, 'color', c(cidx(ii),:))
        end
        Legend{ii} = [conditions{ii} ' ' num2str(mean(cond(ii).dev_angle(15:end)))];
    end
legend(Legend, 'Location', 'eastoutside')
 title('angle')
 
%speed
figure 
    for ii = 1:length(cond)
        hold on
        if ii == 3 || ii == 4 || ii == 7|| ii == 8
            plot(cond(ii).dev_speed,'-.o', 'LineWidth',2, 'color', c(cidx(ii),:))
        else
            plot(cond(ii).dev_speed,'-o', 'LineWidth',2, 'color', c(cidx(ii),:))
        end
        Legend{ii} = [conditions{ii} ' ' num2str(mean(cond(ii).dev_speed(15:end)))];
    end
legend('Location', 'eastoutside')
 title('speed')
 
 %deviation
origin = zeros(1,25);
figure 
    for ii = 1:length(cond)
        hold on
        if ii == 3 || ii == 4 || ii == 7|| ii == 8
            quiver(origin, origin, cond(ii).dev_history(1,:),-cond(ii).dev_history(2,:),'LineWidth',2,'Color', c(cidx(ii),:))
        else
            quiver(origin, origin, cond(ii).dev_history(1,:),-cond(ii).dev_history(2,:),'LineWidth',2,'Color', c(cidx(ii),:))
        end
        
        Legend{ii} = [conditions{ii} ' ' num2str(mean(cond(ii).dev_speed(15:end)))];    
    end
legend('Location', 'eastoutside')

if isfield(cond, 'dev_original')
    for ii = 1:length(cond)
        hold on
        quiver(0, 0, cond(ii).dev_original(1), -cond(ii).dev_original(2), 'Color', 'k')
    end
end
title('deviation')
axis equal

%%
c = flipud(jet(length(cond)));
% conditions = {'down off','down off', 'down on', 'down on'}
conditions = {'up off', 'up off', 'up on','up on'};
% conditions = {'right off', 'right off', 'right on','right on', 'left off','left off', 'left on', 'left on'};
% conditions = {'left off','left off', 'left on', 'left on'};
% conditions = {'right off', 'right off', 'right on','right on'};
idx = [1 2 3 4];
% idx = [5 6 7 8];
% %depth
% Legend = cell(length(cond)/2, 1);
% figure 
%     for ii = 1:length(cond)/2
%         hold on
%         if ii == 3 || ii == 4 || ii == 7|| ii == 8
%             plot(cond(idx(ii)).dev_depth_history,'-o','LineWidth',2,'color', c(ii+4,:))
%         else
%             plot(cond(idx(ii)).dev_depth_history, '-o','LineWidth',2,'color', c(ii+1,:))
%         end
%         Legend{ii} = [conditions{ii} ' ' num2str(mean(cond(idx(ii)).dev_depth_history(15:end)))];
%     end
% legend(Legend, 'Location', 'eastoutside','FontSize',14)
% title('depth','FontSize',14)
% xlabel('trial')
% ylabel('Depth (m)')
% set(gca,'FontSize',12)
 
%angle
figure 
    for ii = 1:length(cond)/2
        hold on
        if ii == 3 || ii == 4 || ii == 7|| ii == 8
            plot(cond(idx(ii)).dev_angle_history, '-o', 'LineWidth',2,'color', c(ii+4,:))
        else
            plot(cond(idx(ii)).dev_angle_history, '-o', 'LineWidth',2,'color', c(ii+1,:))
        end
        Legend{ii} = [conditions{ii} ' ' num2str(mean(cond(idx(ii)).dev_angle_history(end-5:end)))];
    end
legend('Location', 'eastoutside', 'FontSize',14)
 title('angle','FontSize',14)
 xlabel('trial')
 set(gca,'FontSize',12)
 
%speed
figure 
    for ii = 1:length(cond)/2
        hold on
        if ii == 3 || ii == 4 || ii == 7|| ii == 8
            plot(cond(idx(ii)).dev_speed_history, '-o','LineWidth',2,'color', c(ii+4,:))
        else
            plot(cond(idx(ii)).dev_speed_history, '-o', 'LineWidth',2,'color', c(ii+1,:))
        end
        Legend{ii} = [conditions{ii} ' ' num2str(mean(cond(idx(ii)).dev_speed_history(end-5:end)))];
    end
legend('Location', 'eastoutside','FontSize',14)
xlabel('trial')
title('speed','FontSize',14)
set(gca,'FontSize',12)
 
 
%deviation
origin = zeros(1,length(cond(1).dev_history));
figure 
    for ii = 1:length(cond)/2
        hold on
        if ii == 3 || ii == 4 || ii == 7|| ii == 8
            quiver(origin, origin, cond(idx(ii)).dev_history(1,:),-cond(idx(ii)).dev_history(2,:),'LineWidth',2,'Color', c(ii+4,:))
        else
            quiver(origin, origin, cond(idx(ii)).dev_history(1,:),-cond(idx(ii)).dev_history(2,:), 'LineWidth',2,'Color', c(ii+1,:))
        end
        
        Legend{ii} = [conditions{ii}];
    end
% legend(Legend, 'Location', 'eastoutside','FontSize',14)

if isfield(cond, 'original')
    for ii = 1:length(cond)/2
        hold on
        quiver(0, 0, cond(idx(ii)).original(1), -cond(idx(ii)).original(2),'LineWidth',2, 'Color', 'k')
    end
end
title('deviations','FontSize',14)
set(gca,'FontSize',12)
axis equal
% xlim([0 2.5]), ylim([0 2.5])
%%


% -- get colors
c = jet(length(cond));
c = c*.95;

% --- condition names
conditions = {'down off','down off', 'down on', 'down on','up off', 'up off', 'up on','up on'};
% conditions = {'right off', 'right off', 'right on','right on', 'left off','left off', 'left on', 'left on'};

% % % load data into data struct
% % cd './data'
% % initials  = 'hl_pilot';                          % Subject initials used to save file

idx = [16 17];
for ii = 1:length(idx)
    session = idx(ii);
    fname = sprintf ( '%s%i.mat', initials, session);
    load(fname)
    data{ii} = cond;
end

figure
for d = 1:length(data)
    %depth
%     Legend = cell(length(data{1}), 1); 
    for ii = 1:length(data{1})
        hold on
        if ii == 3 || ii == 4 || ii == 7|| ii == 8
            plot(data{d}(ii).dev_depth_history, '-.', 'color', c(ii,:))
        else
            plot(data{d}(ii).dev_depth_history, '-o', 'color', c(ii,:))
        end
%         Legend{ii} = [data{d}(ii) ' ' num2str(mean(data{d}(ii).dev_depth_history(15:end)))];
    end
%     legend(Legend, 'Location', 'eastoutside')
    title('depth')
end

%speed
figure 
    for ii = 1:length(cond)
        hold on
        if ii == 3 || ii == 4 || ii == 7|| ii == 8
            plot(cond(ii).dev_speed, '-.', 'color', c(ii,:))
        else
            plot(cond(ii).dev_speed, '-o', 'color', c(ii,:))
        end
        Legend{ii} = [conditions{ii} ' ' num2str(mean(cond(ii).dev_speed(15:end)))];
    end
legend(Legend, 'Location', 'eastoutside')
 title('speed')
 
 %% Aug 2020
 % -- get colors
c = jet(length(cond));
c = c*.95;

% --- condition names
conditions = {'down off','down off', 'down on', 'down on','up off', 'up off', 'up on','up on'};
% conditions = {'right off', 'right off', 'right on','right on', 'left off','left off', 'left on', 'left on'};

%angle
figure 
    for ii = 1:length(cond)
        hold on
        if cond(ii).dev_index
            plot(cond(ii).dev_angle_history, '-o', 'color', c(ii,:), 'LineWidth',2, 'DisplayName', 'on constraint')
        else
            plot(cond(ii).dev_angle_history, '--', 'color', c(ii,:), 'DisplayName', 'off constraint')
        end
%         Legend{ii} = [conditions{ii} ' ' num2str(mean(cond(ii).dev_angle_history(end-5:end)))];
    end
legend('Location', 'eastoutside')
 title('angle')
 

 
 %deviation
origin = zeros(1,length(cond(1).dev_history));
figure 
    for ii = 1:length(cond)
        hold on
        if cond(ii).dev_index
            quiver(origin, origin, cond(ii).dev_history(1,:),-cond(ii).dev_history(2,:), 'Linestyle','-.', 'Color', c(ii,:), 'AutoScaleFactor',1)
            hold on
            plot(cond(ii).velocity_range(1,:), -cond(ii).velocity_range(2,:),'color',[.5 0 0], 'DisplayName', 'on constraint')
        else
            quiver(origin, origin, cond(ii).dev_history(1,:),-cond(ii).dev_history(2,:), 'Color', c(ii,:), 'AutoScaleFactor',1)
            hold on
            plot(cond(ii).velocity_range(1,:), -cond(ii).velocity_range(2,:),'color',[0 0 .5], 'DisplayName', 'off constraint')
        end
        
%         Legend{ii} = [conditions{ii} ' ' num2str(mean(cond(ii).dev_speed_history(15:end)))];
    end
legend('Location', 'eastoutside')

for ii = 1:length(cond)
    hold on
    quiver(0, 0, cond(ii).original(1), -cond(ii).original(2), 'Color', 'k','AutoScaleFactor',1)
end

title('deviation')
axis equal
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
% hold on, plot(cond(1).velocity_range(1,:), -cond(1).velocity_range(2,:),'color',[.5 0 0], 'DisplayName', 'off constraint range')
% hold on, plot(cond(2).velocity_range(1,:), -cond(2).velocity_range(2,:), 'color',[0 0 .5], 'DisplayName', 'on constraint range')
% 
% axis equal 
% xlim([-5,5]), ylim([-5,5]) 
% legend  
%%
%deviations
figure %(up)
for ii = 1:length(cond(1).dev_history)
    hold on
    quiver(0,0,(cond(1).dev_history(1,ii)), -(cond(1).dev_history(2,ii)), 'color', 'r', 'AutoScaleFactor',1)
    hold on
    quiver(0,0,(cond(2).dev_history(1,ii)), -(cond(2).dev_history(2,ii)), 'color', 'b','AutoScaleFactor',1)
    
end
hold on
quiver(0,0,(cond(1).original(1)), -(cond(1).original(2)), 'color', 'k','AutoScaleFactor',1)
axis equal 

figure %(down)
for ii = 1:length(cond(3).dev_history)
    hold on
    quiver(0,0,abs(cond(3).dev_history(1,ii)), -(cond(3).dev_history(2,ii)), 'color', 'r', 'AutoScaleFactor',1)
    hold on
    quiver(0,0,abs(cond(4).dev_history(1,ii)), -(cond(4).dev_history(2,ii)), 'color', 'b','AutoScaleFactor',1)
end
hold on
quiver(0,0,abs(cond(3).original(1)), -(cond(3).original(2)), 'color', 'k','AutoScaleFactor',1)
axis equal

figure
for ii = 1:length(cond(1).dev_history)
    hold on
    quiver(0,0,(cond(1).dev_history(1,ii)), -(cond(1).dev_history(2,ii)), 'color', 'r', 'AutoScaleFactor',1)
    hold on
    quiver(0,0,(cond(2).dev_history(1,ii)), -(cond(2).dev_history(2,ii)), 'color', 'r','AutoScaleFactor',1)
    hold on
    quiver(0,0,(cond(3).dev_history(1,ii)), -(cond(3).dev_history(2,ii)), 'color', 'b', 'AutoScaleFactor',1)
    hold on
    quiver(0,0,(cond(4).dev_history(1,ii)), -(cond(4).dev_history(2,ii)), 'color', 'b','AutoScaleFactor',1)
end
hold on
quiver(0,0,abs(cond(1).original(1)), -(cond(1).original(2)), 'color', 'k','AutoScaleFactor',1)
axis equal


%angle
figure
for ii = 1:length(cond(1).dev_history)

    hold on
    plot(cond(1).dev_angle_history, '-o','LineWidth',2,'color', 'r')
    hold on
    plot(cond(2).dev_angle_history, '-o','LineWidth',2,'color', 'b')
    hold on
    plot(cond(3).dev_angle_history, '-o','LineWidth',2,'color', 'r')
    hold on
    plot(cond(4).dev_angle_history, '-o','LineWidth',2,'color', 'b')
    
end

%speed
figure
for ii = 1:length(cond(1).dev_history)
    hold on
    plot(cond(1).dev_speed_history, '-o','LineWidth',2,'color', 'r')
    hold on
    plot(cond(2).dev_speed_history, '-o','LineWidth',2,'color', 'b')
    hold on
    plot(cond(3).dev_speed_history, '-o','LineWidth',2,'color', 'r')
    hold on
    plot(cond(4).dev_speed_history, '-o','LineWidth',2,'color', 'b')
end


%% Figures for committee meeting 1
% load hl_pilot12, 16, 14
%deviation

c = flipud(jet(length(cond)));
c1 = brighten(c, .5);
idx = [1 2 3 4];
origin = zeros(1,25);
figure 
    for ii = 1:length(cond)
        hold on
        if  ii == 1 
            scatter(abs(cond(idx(ii)).dev_history(1,:)), -cond(idx(ii)).dev_history(2,:), 50, c1(1,:), 'filled', 'MarkerFaceAlpha', .3)
%             quiver(origin, origin, cond(idx(ii)).dev_history(1,:),-cond(idx(ii)).dev_history(2,:),'LineWidth',2,'Color', c1(7,:), 'ShowArrowHead', 'off', 'AutoScaleFactor',1)

        elseif ii == 2 
            scatter(abs(cond(idx(ii)).dev_history(1,:)), -cond(idx(ii)).dev_history(2,:), 50,  c1(4,:), 'filled', 'MarkerFaceAlpha', .3)
%             quiver(origin, origin, cond(idx(ii)).dev_history(1,:),-cond(idx(ii)).dev_history(2,:),'LineWidth',2,'Color', c1(3,:), 'ShowArrowHead', 'off', 'AutoScaleFactor',1)

        end
        
%         Legend{ii} = [conditions{ii}];
    end
% legend(Legend, 'Location', 'eastoutside','FontSize',14)
axis equal
% xlim([0 2.5]), ylim([0 2.5])
set(gca,'FontSize',14)
title('deviations')


mean1 = mean(abs(cond(idx(2)).dev_history(:,15:end)), 2);
mean2 = mean(abs(cond(idx(2)).dev_history(:,15:end)), 2);
onMean = (mean1 + mean2)./2;
hold on, scatter(onMean(1), -onMean(2), 100, c(4,:), 'filled')
hold on, quiver(0, 0, onMean(1), -onMean(2),'LineWidth',2,'Color', c(4,:), 'ShowArrowHead', 'off', 'AutoScaleFactor',1)

mean1 = mean(abs(cond(idx(1)).dev_history(:,15:end)), 2);
mean2 = mean(abs(cond(idx(1)).dev_history(:,15:end)), 2);
onMean = (mean1 + mean2)./2;
hold on, scatter(onMean(1), -onMean(2), 100,  c(1,:), 'filled')
hold on, quiver(0, 0, onMean(1), -onMean(2),'LineWidth',2,'Color', c(1,:), 'ShowArrowHead', 'off', 'AutoScaleFactor',1)
% 
% axis equal
% xlim([0 5]), ylim([0 5])
% yticks([1 2 3 4 5])
% xticks([1 2 3 4 5])
% xlabel('vx deg/s')
% ylabel('vy deg/s')
% set(gca,'FontSize',14)
%% September 2020
 % -- get colors
c = jet(length(cond));
c = c*.95;

figure %staircase
for ii = 1:length(cond)
    hold on
    if cond(ii).dev_index %if on the constraint
        plot(-cond(ii).step_history, '-o', 'color', c(ii,:), 'LineWidth',2, 'DisplayName', 'on constraint')
    else
        plot(-cond(ii).step_history, '--', 'color', c(ii,:), 'LineWidth', 2, 'DisplayName', 'off constraint')
    end
end
legend
title('staircase steps')

figure %angle
for ii = 1:length(cond)/2
    hold on
    if cond(ii).dev_index %if on the constraint
        plot(cond(ii).dev_angle_history, '-o', 'color', c(ii,:), 'LineWidth',2, 'DisplayName', 'on constraint')
    else
        plot(cond(ii).dev_angle_history, '--', 'color', c(ii,:), 'LineWidth', 2, 'DisplayName', 'off constraint')
    end
end
legend
title('angle steps')

 %deviation
origin = zeros(1,length(cond(1).dev_history));

tested_size = 20;
avg_size = 50;
width = 2; %line width range of velocities and original velocity
figure 
    for ii = 1:length(cond)
        hold on
        if cond(ii).dev_index
            %plot deviations tested
            scatter( cond(ii).dev_history(1,:),-cond(ii).dev_history(2,:),tested_size,c(ii,:), 'filled','DisplayName', 'on constraint')
            
            %plot average of last 5 trials
            avg = mean(cond(ii).dev_history(:,end-4:end),2);
            hold on, plot([0,avg(1)], -[0, avg(2)], 'color',c(ii,:) *.75, 'LineWidth', width-1, 'DisplayName', 'avg last 5 trials')
            
            %calculate distance average to edge of possible velocity range
            %based on depth map
            dist2range = norm(avg - cond(ii).velocity_range(:,2));
            hold on, scatter(avg(1), -avg(2), avg_size, c(ii,:)*.75, 'filled','DisplayName', sprintf('distance: %.2f', dist2range))
            
            %plot velocity range
            hold on
            plot(cond(ii).velocity_range(1,:), -cond(ii).velocity_range(2,:),'color',[1 1 1]*.5, 'LineWidth', width, 'DisplayName', ...
                sprintf('range velocities, translate = %.2f %.2f %.2f', cond(ii).translate))
            

        else
            scatter(cond(ii).dev_history(1,:),-cond(ii).dev_history(2,:), tested_size,c(ii,:), 'filled', 'DisplayName', 'off constraint')
            avg = mean(cond(ii).dev_history(:,end-4:end),2);
            hold on, plot([0,avg(1)], -[0, avg(2)], '--','color',c(ii,:) *.75, 'LineWidth', width-1, 'DisplayName', 'avg last 5 trials')
            dist2range = norm(avg - cond(ii).steps(:,end)); 
            hold on, scatter(avg(1), -avg(2), avg_size, c(ii,:)*.75, 'filled','DisplayName', sprintf('distance: %.2f', dist2range))
            hold on
%             plot(cond(ii).velocity_range(1,:), cond(ii).velocity_range(2,:),'--', 'color',[1 1 1]*.5, 'LineWidth', width,'DisplayName', ...
%                 sprintf('range velocities, translate = %.2f %.2f %.2f', cond(ii).translate))
%             
        end
        
    hold on
    quiver(0, 0, cond(ii).steps(1,end), -cond(ii).steps(2,end), 'Color', 'k','LineWidth', width,'AutoScaleFactor',1, 'DisplayName', '')
    end
legend('Location', 'eastoutside')
title('deviations')

axis equal
% 
%%
myFolder = '/ser/1.1/p1/bonnen/Documents/MATLAB/phaseExperiment/phaseExperiment/data/';
% myFolder = '/Volumes/GoogleDrive/My Drive/opticflow/objectDetection/OpticFlow/phaseExperiment/data/';
filePattern = fullfile(myFolder, '2020-Sep-14*');
matFiles = dir(filePattern);
 
for k = 1:length(matFiles)
    baseFileName = matFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    matdata(k) = load(fullFileName);
end

distances(7).angle = cond(1).angle_btwn_constraints;
distances(8).angle = cond(5).angle_btwn_constraints;

for ii = 3:length(matdata)
    for jj = 1:length(matdata(3).cond)
        if matdata(ii).cond(jj).dev_index % if on the constraint
            avg = mean(matdata(ii).cond(jj).dev_history(:,end-4:end),2);
            diff = avg - matdata(ii).cond(jj).velocity_range(:,1);
            dirConstraint = matdata(ii).cond(jj).velocity_range(:,1)-matdata(ii).cond(jj).velocity_range(:,2);
            dist2range = sign(dot(dirConstraint,diff))* norm(diff); %signed distace to edge of range
            if round(distances(7).angle) == round(matdata(ii).cond(jj).angle_btwn_constraints) 
                distances(7).on = [distances(7).on, dist2range];
            else
                distances(8).on = [distances(8).on, dist2range];
            end
                        % calculate constraint line segment
            range_edge_near = matdata(ii).cond(jj).velocity_range(:,1);
            range_edge_far = matdata(ii).cond(jj).velocity_range(:,2);
            range_vec = range_edge_near-range_edge_far;
            x1 = range_edge_near(1);
            x2 = range_edge_far(1);
            y1 = range_edge_near(2);
            y2 = range_edge_far(2);
            m1 = (y2-y1)/(x2-x1);
            b1 = range_edge_near(2)-m1*range_edge_near(1);
        else
            avg = mean(matdata(ii).cond(jj).dev_history(:,end-4:end),2);
            % calculate distance to line segment
            % dist2line = abs((y2-y1)*avg(1) - (x2-x1)*avg(2) + x2*y1-y2*x1)/sqrt((y2-y1)^2+(x2-x1)^2);
          
            m2 = -1/m1;
            b2 = avg(2)-(m2)*avg(1);
            intersect = [(b2-b1)/(m1-m2); m1*(b2-b1)/(m1-m2)+b1];
            dist2line = norm(avg - intersect);
            dist2near = norm(intersect - range_edge_near);
            dist2far = norm(intersect - range_edge_far);
            if (dist2near && dist2far) <= norm(range_edge_near - range_edge_far) %if orthogonal projection lies within segment
                if round(distances(7).angle) == round(matdata(ii).cond(jj).angle_btwn_constraints) %first translation
                    distances(7).off = [distances(7).off, dist2range];
                else % second translation
                    distances(8).off = [distances(8).off, dist2range];
                end
            else
                %set as distance to closest edge
                if round(distances(7).angle) == round(matdata(ii).cond(jj).angle_btwn_constraints)
                    distances(7).off = [distances(7).off, min(norm(avg - range_edge_near), norm(avg - range_edge_far))];
                else % second translation
                    distances(8).off = [distances(8).off, min(norm(avg - range_edge_near), norm(avg - range_edge_far))];
                end
            end
        end
    end
end
      
 %% rotating constraint 
% myFolder = '/ser/1.1/p1/bonnen/Documents/MATLAB/phaseExperiment/phaseExperiment/data/';
myFolder = '/Volumes/GoogleDrive/My Drive/opticflow/objectDetection/OpticFlow/phaseExperiment/data/archive';
% filePattern = fullfile(myFolder, '2020-Sep-17_HL_ground_plane*'); %includes 133
% filePattern = fullfile(myFolder, '2020-Sep-18_HL_GP_45*');
filePattern = fullfile(myFolder, '2020-Sep-14_HL_*'); %archive

% filePattern = fullfile(myFolder, '2020-Sep-26_HL_noise_gp77*');
% filePattern = fullfile(myFolder, '2020-Sep-27_HL_noise_gp45*');
% filePattern = fullfile(myFolder, '2020-Sep-27_HL_noise_gp77*');

matFiles = dir(filePattern);
 
for k = 1:length(matFiles)
    baseFileName = matFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    matdata(k) = load(fullFileName);
end

angles = [0 45 90 133 135 180 225 270 315];
for ii = 1:length(matdata)+1
    distances(ii).angle = angles(ii);
    distances(ii).on = [];
    distances(ii).off = [];
    
    velocities(ii).angle = angles(ii);
    velocities(ii).on = [];
    velocities(ii).off = [];
end

grad = linspace(1, 1, length(matdata));
figure(1)
for ii = 1:length(matdata)
    for jj = 1:length(matdata(1).cond)
        if matdata(ii).cond(jj).dev_index % if on the constraint
            avg = mean(matdata(ii).cond(jj).dev_history(:,end-4:end),2);
            
            diff = avg - matdata(ii).cond(jj).velocity_range(:,1);
            dirConstraint = matdata(ii).cond(jj).velocity_range(:,1)-matdata(ii).cond(jj).velocity_range(:,2);
            dist2range = sign(dot(dirConstraint,diff))* norm(diff); %signed distace to edge of range

            %save distance
            distances(1).on = [distances(1).on, dist2range];
            velocities(1).on = [velocities(1).on, avg];
                        % calculate constraint line segment
            range_edge_near = matdata(ii).cond(jj).velocity_range(:,1);
            range_edge_far = matdata(ii).cond(jj).velocity_range(:,2);
            range_vec = range_edge_near-range_edge_far;
            x1 = range_edge_near(1);
            x2 = range_edge_far(1);
            y1 = range_edge_near(2);
            y2 = range_edge_far(2);
            m1 = (y2-y1)/(x2-x1);
            b1 = range_edge_near(2)-m1*range_edge_near(1);
           
            hold on, scatter(avg(1), -avg(2), 50,[0 1 1]*grad(ii), 'filled')

        else
            avg = mean(matdata(ii).cond(jj).dev_history(:,end-4:end),2);
            % calculate distance to line segment
            % dist2line = abs((y2-y1)*avg(1) - (x2-x1)*avg(2) + x2*y1-y2*x1)/sqrt((y2-y1)^2+(x2-x1)^2);
          
            m2 = -1/m1;
            b2 = avg(2)-(m2)*avg(1);
            intersect = [(b2-b1)/(m1-m2); m1*(b2-b1)/(m1-m2)+b1];
            dist2line = norm(avg - intersect);
            dist2near = norm(intersect - range_edge_near);
            dist2far = norm(intersect - range_edge_far);
            if (dist2near && dist2far) <= norm(range_edge_near - range_edge_far) %if orthogonal projection lies within segment
                distances(ii+1).off = [distances(ii+1).off, dist2range];
            else
                distances(ii+1).off = [distances(ii+1).off, min(norm(avg - range_edge_near), norm(avg - range_edge_far))];
            end
            velocities(ii+1).off = [velocities(ii+1).off, avg];
            hold on, scatter(avg(1), -avg(2), 50, [1 0 1]*grad(ii), 'filled')
            
        end
    end
end

figure
hold on, plot(matdata(1).cond(1).velocity_range(1,:), -matdata(1).cond(1).velocity_range(2,:),'color',[1 1 1]*.5)
hold on, quiver(0, 0, matdata(1).cond(1).steps(1,end), -matdata(1).cond(1).steps(2,end), 'Color', [.2, .2, .7],'AutoScaleFactor',1)
axis equal

xlim([0, 5])
ylim([-2.5,2.5])
axis square

% figure(2)
% polarplot(0, mean(distances(1).on), 'b-o')
% hold on
% for a = 2:length(angles)
%     hold on
%     polarplot(deg2rad(angles(a)), mean(distances(a).off), 'r-o')
% end
% 
figure
% plot actual means
thresh = zeros(2,7);
for v = 2:length(velocities)
    avg = mean(velocities(v).off,2);
    
    thresh(:,v-1) = avg;
    hold on, scatter(avg(1), -avg(2), 50,[.2 .2 .7]*grad(ii), 'filled')
end
avg = mean(velocities(1).on,2);
thresh = [thresh avg thresh(:,1)];
hold on, scatter(avg(1), -avg(2), 50,[.2 .2 .7]*grad(ii), 'filled')
hold on, plot(matdata(1).cond(1).velocity_range(1,:), -matdata(1).cond(1).velocity_range(2,:),'color',[1 1 1]*.5)
hold on, quiver(0, 0, matdata(1).cond(1).steps(1,end), -matdata(1).cond(1).steps(2,end), 'Color', [.2 .2 .7],'LineWidth', 2,'AutoScaleFactor',1)
hold on, plot(thresh(1,:), -thresh(2,:), 'color', [.2 .2 .7], 'linewidth', 2)
axis equal
xlim([0, 5])
ylim([-2.5,2.5])
axis square

% % flow parsing
% % new_center = matdata(1).cond(1).steps(:,end);
% % 
% % % plot recentered means
% figure(4)
% current_center = matdata(1).cond(1).steps(:,end);
% for v = 2:length(velocities)
%     avg = mean(velocities(v).off,2)-current_center;
%     hold on, scatter(avg(1), -avg(2), 50,[0 1 1]*grad(ii), 'filled')
% end
% avg = mean(velocities(1).on,2)-current_center;
% hold on, scatter(avg(1), -avg(2), 50,[0 1 1]*grad(ii), 'filled')
% new_range = matdata(1).cond(1).velocity_range -current_center;
% hold on, plot(new_range(1,:), -new_range(2,:),'color',[1 1 1]*.2)
% axis equal
% grid on

% get ellipse
clear all_avg
all_avg = mean(velocities(1).on,2);
for ii = 2:length(velocities)
    avg = mean(velocities(ii).off,2);
    all_avg = [all_avg avg];
end
all_avg = all_avg';
center = mean(all_avg);
demeaned_all_avg = all_avg-center;

C=demeaned_all_avg'*demeaned_all_avg/(14);%find the covariance matrix
[V,D]=eig(C);
[d,ind] = sort(diag(D), 'descend');
D = D(ind,ind);
V = V(:,ind);
v1=V(:,1)*(D(1,1));%eigenvectors
v2=V(:,2)*(D(2,2));

% hold on
% plot([center(1), center(1)+v1(1)],-[center(2), center(2)+v1(2)],'b-', 'LineWidth', 2, 'DisplayName', 'PC 1');
% hold on
% plot([center(1), center(1)+v2(1)],-[center(2), center(2)+v2(2)],'r-', 'LineWidth', 2, 'DisplayName', 'PC 1');

theta = linspace(0, 2*pi, 50)';
x = [cos(theta) sin(theta)];
ellipse = 5*[v1 v2]*x';
addmean = ellipse+center';

% figure(10)
% hold on
% plot([center(1), center(1)+1/s(1,1)*v1(1)],-[center(2),center(2)+1/s(1,1)*v1(2)])
% hold on
% plot([center(1), center(1)+1/s(2,2)*v(1,1)],-[center(2),center(2)+1/s(2,2)*v(2,1)])

% plot ellipse and mean
% hold on
% plot(addmean(1,:), -addmean(2,:), 'color', [.2 .2 .7], 'Linewidth', 2)
% hold on, scatter(center(1), -center(2), 50, [.2 .2 .7], 'filled')
% hold on
% scatter(all_avg(:,1), -all_avg(:,2))

 %%
figure
for ii = 1:length(distances)
    hold on, scatter(ones(1,length(distances(ii).on))*distances(ii).angle, distances(ii).on, 'b')
    hold on, plot([distances(ii).angle-5, distances(ii).angle+5], [mean(distances(ii).on), mean(distances(ii).on)], 'b', 'LineWidth', 2)
    hold on, plot([distances(ii).angle, distances(ii).angle], [mean(distances(ii).on)- std(distances(ii).on), mean(distances(ii).on)+std(distances(ii).on)], 'color', [0 0 1], 'LineWidth', 2)
    hold on, scatter(ones(1,length(distances(ii).on))*distances(ii).angle, distances(ii).off,'r')
    hold on, plot([distances(ii).angle-5, distances(ii).angle+5], [mean(distances(ii).off), mean(distances(ii).off)], 'r', 'LineWidth', 2)
    hold on, plot([distances(ii).angle, distances(ii).angle], [mean(distances(ii).off)- std(distances(ii).off), mean(distances(ii).off)+std(distances(ii).off)], 'color',[1 0 0], 'LineWidth', 2)
end


title('distance to constraint segment vs angle between on and off constraint')
ax = gca;
ax.FontSize = 16;
xlim([0 180])
ylabel('distance to constraint segment (velocity space)')
xlabel('angle between on and off constraint (deg)')


%% Nov 2020 --- controls

 % -- get colors
c = jet(length(cond));
c = c*.95;

%plot staircases
angles = [0 45 90 135 180 225 270 315 0 180];
num_staircases = 2;
num_t = 1;


% figure
% for ii = 1:length(cond)/num_staircases
%      subplot(num_t,8,ii)
%      plot(-cond(2*ii-1).step_history, '-o', 'color', c(2*ii-1,:), 'LineWidth',2, 'DisplayName', num2str(cond(2*ii-1).rotation))
%      hold on
%      plot(-cond(2*ii).step_history, '-o', 'color', c(2*ii,:), 'LineWidth',2, 'DisplayName', num2str(cond(2*ii).rotation))
%      title(num2str(cond(2*ii).rotation))
%      ylim([-10 0])
% end

% title('staircase steps')


% calculate averages
avg = zeros(2,length(angles)*num_t);
for ii = 1:length(cond)
    avg(:,ii) = mean(cond(ii).dev_history(:,end-4:end),2);
end       
a = mean(reshape(avg', 2,length(angles)*num_t,2), 1);
means = [a(:,:,1); a(:,:,2)];

figure
scatter(cond(1).steps(1,end), -cond(1).steps(2,end), 50, 'k', 'filled')

% hold on
% plot([means(1,1:length(angles)), means(1,1)], -[means(2,1:length(angles)),means(2,1)], 'o-','color', [.2, .2, .7], 'linewidth', 2)
% hold on
% plot([means(1,length(angles)+1:end), means(1,9)], -[means(2,length(angles)+1:end),means(2,9)], 'o-','color', [.2, .7, .2], 'linewidth', 2)
% hold on
% quiver(0,0,cond(1).steps(1,end), -cond(1).steps(2,end), 'Color', [.5, .5, .7],'LineWidth', 2,'AutoScaleFactor',1)
% hold on
% quiver(0,0,cond(17).steps(1,end), -cond(17).steps(2,end), 'Color', [.5, .7, .5],'LineWidth', 2,'AutoScaleFactor',1)
% xlim([0 3])
% ylim([-3 0])
% axis equal
% 

hold on
scatter(means(1,1:length(angles)), -means(2,1:length(angles)), 50, [.2, .2, .7], 'filled')
% xlim([0 3])
% ylim([-3 0])
axis equal

 allX = [allX; means(1,:)];
 allY = [allY; means(2,:)];

%%
% 
% allX = [allmeans(1, 1:8);allmeans(1,9:16)];
% allY = [allmeans(2, 1:8);allmeans(2,9:16)];
% 
% allX(3:4,:) = -[allX(3:4,1) flip(allX(3:4,2:end),2)];
% allY(3:4,:) = [allY(3:4,1) flip(allY(3:4,2:end),2)];

% allX = [allmeans(1, 1:8);allmeans(1,9:16); -[allmeans(1,17), flip(allmeans(1,18:24))]; -[allmeans(1,25), flip(allmeans(1,26:32))]];
% allY = [allmeans(2, 1:8);allmeans(2,9:16); [allmeans(2,17), flip(allmeans(2,18:24))]; [allmeans(2,25), flip(allmeans(2,26:32))]];
meanX = mean((allX));
meanY = mean(allY);

% flatX = allX(:);
% flatY = allY(:);
% figure
% scatter(flatX(1:16), -flatY(1:16),50,[.2, .2, .7], 'filled')
% hold on
% scatter(flatX(17:end), -flatY(17:end),50,[.2, .7, .2], 'filled')
% 

figure
hold on
plot([meanX(1:8) meanX(1)], -[meanY(1:8) meanY(1)], 'color',[.7, .5, .5], 'linewidth', 2)
hold on, quiver(0, 0, cond(1).steps(1,end), -cond(1).steps(2,end), 'Color', 'k','LineWidth', 2,'AutoScaleFactor',1)
axis equal

figure
hold on
plot([meanX(9:end) meanX(9)], -[meanY(9:end) meanY(9)],'color',[.7, .5, .5], 'linewidth', 2)
hold on, quiver(0, 0, cond(17).steps(1,end), -cond(17).steps(2,end), 'Color', 'k','LineWidth', 2,'AutoScaleFactor',1)
% xlim([0 5])
% ylim([-2.5 2.5])
axis equal

hold on
for c = 1:length(cond)/2
    scatter(cond(c).steps(1,1), -cond(c).steps(2,1), 50, 'k', 'filled')
end

idx = [[1:8];[9:16]];
idx = 1:10
colors = [[.2 .2 .7];[.2 .7 .2]];
figure
for t = 1:num_t
    meanX = means(1,idx(t,:));
    meanY = means(2,idx(t,:));
    center = [mean(meanX), mean(meanY)];
    demeaned_all_avg = [meanX' meanY']-center;

    C=demeaned_all_avg'*demeaned_all_avg/(14);%find the covariance matrix
    [V,D]=eig(C);
    [d,ind] = sort(diag(D), 'descend');
    D = D(ind,ind);
    V = V(:,ind);
    v1=V(:,1)*(D(1,1));%eigenvectors
    v2=V(:,2)*(D(2,2));
    theta = linspace(0, 2*pi, 50)';
    x = [cos(theta) sin(theta)];
    ellipse = 5.5*[v1 v2]*x';
    addmean = ellipse+center';

    hold on
    plot(addmean(1,:), -addmean(2,:), 'color', colors(t,:), 'Linewidth', 2)
    hold on, scatter(center(1), -center(2), 50, colors(t,:), 'filled')
end
%weber fraction for in same direction
% weber = (abs(vecnorm(means) - norm(cond(1).steps(:,end))))/ norm(cond(1).steps(:,end))
% 
% thresh = [meanX;meanY];
% weber = (abs(vecnorm(thresh) -  norm(cond(1).steps(:,end))))/ norm(cond(1).steps(:,end))
%          