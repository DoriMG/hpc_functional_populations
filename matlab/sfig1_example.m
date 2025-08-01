
% Example plot
subplot(1,2,1);
passive = [linspace(0,1,8), ones(1,6)]+(rand(1,14)-0.5)/4;
passive(passive>1) = 1;
passive(passive<0) = 0;
passive_control =  [linspace(0,0.4,7), linspace(0.4,0.4,7)]+(rand(1,14)-0.5)/4;
passive_control(passive_control>1) = 1;
passive_control(passive_control<0) = 0;
plot(passive_control);hold on;
plot(passive);
% legend({'control trials', 'rewarded trials'})
ylabel('Ratio trials licked')
title('Passive sessions')

subplot(1,2,2);
active = [linspace(0,0.1,4), linspace(0.1,0.7,7), ones(1,3)*0.7]+(rand(1,14)-0.5)/4;
active(active>1) = 1;
active(active<0) = 0;
active_control = [linspace(0,0.1,4), linspace(0.1,0.3,7), ones(1,3)*0.3]+(rand(1,14)-0.5)/4;
active_control(active_control>1) = 1;
active_control(active_control<0) = 0;
plot(active_control);hold on;
plot(active);
legend({'control trials', 'rewarded trials'}, 'Location', 'southeast')
ylabel('Ratio trials licked')
title('Active sessions')


datatemp = [passive, passive_control, active, active_control];
ap = [ones(size(passive)), ones(size(passive_control)), ones(size(active))*2, ones(size(active_control))*2];
sessions = [[1:length(passive)], [1:length(passive_control)], [1:length(active)], [1:length(active_control)]];
control = [ones(size(passive)), ones(size(passive_control))*2, ones(size(active)), ones(size(active_control))*2];
data = [datatemp; sessions;control; ap]';
save('R/data/example_behaviour.mat', 'data')