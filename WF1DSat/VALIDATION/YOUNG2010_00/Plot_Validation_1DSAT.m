%PLOT VALIDATION 1DSAT
%% Read inputs on nodes from FLOW1DUNSAT
inputfilename = 'sat_inputs.wfsinp';
[p,f,e]=fileparts(inputfilename);
inputfilename=fullfile(p,f);
satoutputnodes = readtable(strcat(inputfilename,'.outnods'),'FileType','delimitedtext');

% satoutputnodes = readtable('sat_outputs_nodes.csv');
HeadVariation_1layer_Young2010 = readtable('DATA_YOUNG_2010\Head_Variation_100.csv');
Watertables_1layer_Young2010 = readtable('DATA_YOUNG_2010\Watertables.csv');
Outflow_1layer_Young2010 = readtable('DATA_YOUNG_2010\Outflow.csv');


%% Plot outflow on time
x0 = 0.0;
x1 = 1.0;
x2 = 2.0;
xout = 2.25;

close all
figure1 = figure('Renderer', 'painters', 'Position', [0 0 800 300],'Color',[1 1 1]);

axes1 = axes('Parent',figure1);
hold(axes1,'on');

p1 = plot(Outflow_1layer_Young2010.t,Outflow_1layer_Young2010.Outflow,'k');
hold on
s1 = scatter(satoutputnodes.t(satoutputnodes.x==xout),satoutputnodes.qhor(satoutputnodes.x==xout),'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1);
s1.SizeData = 10.0;
%%Plot format
xlabel('time (days)') 
ylabel('outflow (m3/d)') 
legend([p1,s1],{'Young,2010','FLOW1DSAT'},'Location','northeast')
set(axes1,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')
hold off

%% Plot waterhead in time
figure2 = figure('Renderer', 'painters', 'Position', [0 0 800 300],'Color',[1 1 1]);
axes2 = axes('Parent',figure2);
hold(axes2,'on');
p21 = plot(HeadVariation_1layer_Young2010.t,HeadVariation_1layer_Young2010.Head0,'k');
p22 = plot(HeadVariation_1layer_Young2010.t,HeadVariation_1layer_Young2010.Head1,'k');
p23 = plot(HeadVariation_1layer_Young2010.t,HeadVariation_1layer_Young2010.Head2,'k');
p24 = plot(HeadVariation_1layer_Young2010.t,HeadVariation_1layer_Young2010.HeadBoundary,'k');

s21 = scatter(satoutputnodes.t(satoutputnodes.x==x0),satoutputnodes.head(satoutputnodes.x==x0),'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1);
s22 = scatter(satoutputnodes.t(satoutputnodes.x==x1),satoutputnodes.head(satoutputnodes.x==x1),'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1);
s23 = scatter(satoutputnodes.t(satoutputnodes.x==x2),satoutputnodes.head(satoutputnodes.x==x2),'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1);
s24 = scatter(satoutputnodes.t(satoutputnodes.x==xout),satoutputnodes.head(satoutputnodes.x==xout),'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1);

%%Plot format
xlabel('time (days)') 
ylabel('water table height (m)') 
p21.LineWidth = 1.5;
p22.LineWidth = 1.5;
p23.LineWidth = 1.5;
p24.LineWidth = 1.5;

s21.SizeData = 7.0;
s22.SizeData = 7.0;
s23.SizeData = 7.0;
s24.SizeData = 7.0;
legend([p21,s21],{'Young,2010','FLOW1DSAT'},'Location','northeast')
axis(axes2,[0 0.12 0 0.5])
set(axes2,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')

%% Plot watertables in time
figure3 = figure('Renderer', 'painters', 'Position', [0 0 800 300],'Color',[1 1 1]);
axes3 = axes('Parent',figure3);
hold(axes3,'on');
p31 = plot(Watertables_1layer_Young2010.x,Watertables_1layer_Young2010.h_recharge,'k');
p32 = plot(Watertables_1layer_Young2010.x,Watertables_1layer_Young2010.h_0_2_days,'k');

s31 = scatter(satoutputnodes.x(satoutputnodes.t==0.042),satoutputnodes.head(satoutputnodes.t==0.042),'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1);
s32 = scatter(satoutputnodes.x(satoutputnodes.t==0.2),satoutputnodes.head(satoutputnodes.t==0.2),'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1);
hold on

%%Plot format
xlabel('distance (m)') 
ylabel('water table height (m)') 
p31.LineWidth = 1.5;
p32.LineWidth = 1.5;
s31.SizeData = 7.0;
s32.SizeData = 7.0;
legend([p31,s31],{'Young,2010','FLOW1DSAT'},'Location','northeast')
axis(axes3,[0 2.25 0 0.5])
set(axes3,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')
