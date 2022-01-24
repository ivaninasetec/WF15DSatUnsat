%PLOT VALIDATION 1DSAT
%% Read inputs on nodes from FLOW1DUNSAT
inputfilename = 'sat_inputs.wfsinp';
[p,f,e]=fileparts(inputfilename);
inputfilename=fullfile(p,f);
satoutputnodes = readtable(strcat(inputfilename,'.outnods'),'FileType','delimitedtext');
%%satoutputnodes = readtable('sat_outputs_nodes.csv');

%% Plot outflow on time
x0 = 0.0;
x1 = 1.0;
x2 = 2.0;
xout = 2.25;

close all
%% Plot watertables in time
figure3 = figure('Renderer', 'painters', 'Position', [0 0 600 200],'Color',[1 1 1]);
axes3 = axes('Parent',figure3);
hold(axes3,'on');
l0 = plot([0,2.25],[0,0],'k');
l1 = plot([0,2.25],[0.13,0.13],'k');
l2 = plot([0,2.25],[0.3,0.3],'k');
l0.LineWidth = 1.5;
l0.LineWidth = 1.5;

pp11 = plot(satoutputnodes.x(satoutputnodes.t==0.001),satoutputnodes.head(satoutputnodes.t==0.001),'b');
pp12 = plot(satoutputnodes.x(satoutputnodes.t==0.002),satoutputnodes.head(satoutputnodes.t==0.002),'b');
pp13 = plot(satoutputnodes.x(satoutputnodes.t==0.003),satoutputnodes.head(satoutputnodes.t==0.003),'b');
pp14 = plot(satoutputnodes.x(satoutputnodes.t==0.004),satoutputnodes.head(satoutputnodes.t==0.004),'b');
pp15 = plot(satoutputnodes.x(satoutputnodes.t==0.005),satoutputnodes.head(satoutputnodes.t==0.005),'b');
pp16 = plot(satoutputnodes.x(satoutputnodes.t==0.01),satoutputnodes.head(satoutputnodes.t==0.01),'b');
pp17 = plot(satoutputnodes.x(satoutputnodes.t==0.015),satoutputnodes.head(satoutputnodes.t==0.015),'b');
pp18 = plot(satoutputnodes.x(satoutputnodes.t==0.02),satoutputnodes.head(satoutputnodes.t==0.02),'b');

pp19 = plot(satoutputnodes.x(satoutputnodes.t==0.05),satoutputnodes.head(satoutputnodes.t==0.05),'--b');
pp20 = plot(satoutputnodes.x(satoutputnodes.t==0.1),satoutputnodes.head(satoutputnodes.t==0.1),'--b');
pp21 = plot(satoutputnodes.x(satoutputnodes.t==0.15),satoutputnodes.head(satoutputnodes.t==0.15),'--b');
pp22 = plot(satoutputnodes.x(satoutputnodes.t==0.2),satoutputnodes.head(satoutputnodes.t==0.2),'--b');

% s31 = scatter(satoutputnodes.x(satoutputnodes.t==0.001),satoutputnodes.head(satoutputnodes.t==0.001),'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1);
% s31 = scatter(satoutputnodes.x(satoutputnodes.t==0.005),satoutputnodes.head(satoutputnodes.t==0.005),'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1);
% s32 = scatter(satoutputnodes.x(satoutputnodes.t==0.01),satoutputnodes.head(satoutputnodes.t==0.01),'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1);
% s32 = scatter(satoutputnodes.x(satoutputnodes.t==0.02),satoutputnodes.head(satoutputnodes.t==0.02),'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1);
% s31.SizeData = 7.0;
% s32.SizeData = 7.0;

%Plot format
xlabel('distance (m)','Fontsize',8) 
ylabel('water table height (m)')
legend([pp11,pp19],{'Recharge: 0.00d-0.02d','Discharge: 0.02d-0.2d'},'Location','northeast','Fontsize',8)
axis(axes3,[0 2.25 0 0.3])
set(axes3,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')
set(gca,'FontSize',8)
%% Plot watertables in time (on 0.12h and 0.5hr as in Roushton and Youngs)
figure3 = figure('Renderer', 'painters', 'Position', [0 0 250 250],'Color',[1 1 1]);
axes3 = axes('Parent',figure3);
hold(axes3,'on');
l0 = plot([0,2.25],[0,0],'k');
l1 = plot([0,2.25],[0.13,0.13],'k');
l2 = plot([0,2.25],[0.3,0.3],'k');
l0.LineWidth = 1.5;
l0.LineWidth = 1.5;

pp15 = plot(satoutputnodes.x(satoutputnodes.t==0.005),satoutputnodes.head(satoutputnodes.t==0.005),'--b');
pp16 = plot(satoutputnodes.x(satoutputnodes.t==0.01),satoutputnodes.head(satoutputnodes.t==0.021),'b');

% s31 = scatter(satoutputnodes.x(satoutputnodes.t==0.001),satoutputnodes.head(satoutputnodes.t==0.001),'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1);
% s31 = scatter(satoutputnodes.x(satoutputnodes.t==0.005),satoutputnodes.head(satoutputnodes.t==0.005),'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1);
% s32 = scatter(satoutputnodes.x(satoutputnodes.t==0.01),satoutputnodes.head(satoutputnodes.t==0.01),'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1);
% s32 = scatter(satoutputnodes.x(satoutputnodes.t==0.02),satoutputnodes.head(satoutputnodes.t==0.02),'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1);
% s31.SizeData = 7.0;
% s32.SizeData = 7.0;

%Plot format
xlabel('distance (m)') 
ylabel('water table height (m)')
legend([pp16,pp15],{'Watertable at 0.50h','Watertable at 0.12h'},'Location','southwest')
axis(axes3,[1.875 2.25 0 0.16])
set(axes3,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')
set(gca,'FontSize',8)