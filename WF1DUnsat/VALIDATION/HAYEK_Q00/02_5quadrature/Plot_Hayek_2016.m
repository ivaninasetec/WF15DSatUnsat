%SCRIPT FOR THE VALIDATION OF FLOW1DUNSAT

%% Read inputs on nodes from FLOW1DUNSAT
inputfilename = 'unsat_inputs.wfuinp';
[p,f,e]=fileparts(inputfilename);
inputfilename=fullfile(p,f);
unsatoutputnodes = readtable(strcat(inputfilename,'.outnodu'),'FileType','delimitedtext');
% unsatoutputnodes = readtable('unsat_outputs_nodes.csv');

%% Construct object for analytical Hayek, 2016 expressions
% a=0.1cm^-1 n=3.5, k=1cm/h, thres=0.06 thsat=0.4
hayekobj = hayek(0.1,3.5,1.0,0.06,0.4);

%%
% Plot from z0=0 to zf=2m at time 0
close all
% Create figure
figure1 = figure('Color',[1 1 1]);
set(figure1, 'PaperUnits', 'centimeters');
set(figure1, 'PaperSize', [4 4]);
% set(figure1, 'PaperPositionMode', 'manual');
% set(figure1, 'PaperPosition', [0 0 4 4]);
figure1.Position = [100 100 400 300];
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

fig1 = hayekobj.plot_th_in_z_at_t([0:-0.01:-200],3);
hold on
hayekobj.plot_th_in_z_at_t([0:-0.01:-200],6);
hayekobj.plot_th_in_z_at_t([0:-0.01:-200],12);
hayekobj.plot_th_in_z_at_t([0:-0.01:-200],24);
hayekobj.plot_th_in_z_at_t([0:-0.01:-200],36);
hayekobj.plot_th_in_z_at_t([0:-0.01:-200],48);
% s1 = scatter(unsatoutputnodes.th(unsatoutputnodes.t==3),unsatoutputnodes.x(unsatoutputnodes.t==3),'k');
% s2 = scatter(unsatoutputnodes.th(unsatoutputnodes.t==6),unsatoutputnodes.x(unsatoutputnodes.t==6),'filled','k');
% s3 = scatter(unsatoutputnodes.th(unsatoutputnodes.t==12),unsatoutputnodes.x(unsatoutputnodes.t==12),'filled','k');
% s4 = scatter(unsatoutputnodes.th(unsatoutputnodes.t==24),unsatoutputnodes.x(unsatoutputnodes.t==24),'filled','k');
% s5 = scatter(unsatoutputnodes.th(unsatoutputnodes.t==36),unsatoutputnodes.x(unsatoutputnodes.t==36),'filled','k');
% s6 = scatter(unsatoutputnodes.th(unsatoutputnodes.t==48),unsatoutputnodes.x(unsatoutputnodes.t==48),'filled','k');
% s1.SizeData = 10.0;
% s2.SizeData = 10.0;
% s3.SizeData = 10.0;
% s4.SizeData = 10.0;
% s5.SizeData = 10.0;
% s6.SizeData = 10.0;
s1=plot(unsatoutputnodes.th(unsatoutputnodes.t==3),unsatoutputnodes.x(unsatoutputnodes.t==3),'--ok','MarkerSize',3);
plot(unsatoutputnodes.th(unsatoutputnodes.t==6),unsatoutputnodes.x(unsatoutputnodes.t==6),'--ok','MarkerSize',3);
plot(unsatoutputnodes.th(unsatoutputnodes.t==12),unsatoutputnodes.x(unsatoutputnodes.t==12),'--ok','MarkerSize',3);
plot(unsatoutputnodes.th(unsatoutputnodes.t==24),unsatoutputnodes.x(unsatoutputnodes.t==24),'--ok','MarkerSize',3);
plot(unsatoutputnodes.th(unsatoutputnodes.t==36),unsatoutputnodes.x(unsatoutputnodes.t==36),'--ok','MarkerSize',3);
plot(unsatoutputnodes.th(unsatoutputnodes.t==48),unsatoutputnodes.x(unsatoutputnodes.t==48),'--ok','MarkerSize',3);

%%Plot format
xlabel('water content (cm3/cm3)') 
ylabel('depht (cm)') 
legend([fig1,s1],{'Hayek,2016 (3,6,12,24,36,48h)','WF1DUNSAT 3,6,12,24,36,48h)'},'Location','southeast','FontSize',8)
set(axes1,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')

hold off
%% RMSE
time=3;
RMSE3 =sqrt(mean((hayekobj.th_zvec_tsca(unsatoutputnodes.x(unsatoutputnodes.t==time),time)-unsatoutputnodes.th(unsatoutputnodes.t==time)).^2));
time=6;
RMSE6 =sqrt(mean((hayekobj.th_zvec_tsca(unsatoutputnodes.x(unsatoutputnodes.t==time),time)-unsatoutputnodes.th(unsatoutputnodes.t==time)).^2));
time=12;
RMSE12 =sqrt(mean((hayekobj.th_zvec_tsca(unsatoutputnodes.x(unsatoutputnodes.t==time),time)-unsatoutputnodes.th(unsatoutputnodes.t==time)).^2));
time=24;
RMSE24 =sqrt(mean((hayekobj.th_zvec_tsca(unsatoutputnodes.x(unsatoutputnodes.t==time),time)-unsatoutputnodes.th(unsatoutputnodes.t==time)).^2));
time=48;
RMSE48 =sqrt(mean((hayekobj.th_zvec_tsca(unsatoutputnodes.x(unsatoutputnodes.t==time),time)-unsatoutputnodes.th(unsatoutputnodes.t==time)).^2));
RMSE=mean([RMSE3,RMSE6,RMSE12,RMSE24,RMSE48]);
% Create textbox
annotation('textbox',...
    [0.635 0.3 0.24 0.06],...
    'String',{'RMSE='+compose("%8.4f",RMSE)},...
    'FitBoxToText','off','FontSize',8);
