%VALIDATION AGAINST HEYNS

close all
inputfilename = 'Heyns_validation.wfuinp';
[p,f,e]=fileparts(inputfilename);
inputfilename=fullfile(p,f);
unsatoutputelements = readtable(strcat(inputfilename,'.outelmu.csv'),'FileType','delimitedtext');
satoutputelements = readtable(strcat(inputfilename,'.outelms.csv'),'FileType','delimitedtext');
modelconstraints = readtable(strcat(inputfilename,'.outcons.csv'),'FileType','delimitedtext');

% unsatoutputelements = readtable('unsat_outputs_elements.csv');
% satoutputelements   = readtable('sat_outputs_elements.csv');
% modelconstraints    = readtable('model_constraints.csv');
heynsp3 = readtable('HEYNS_DATA\P3.csv');
heynsp5 = readtable('HEYNS_DATA\P5.csv');
HeynsQ = readtable('HEYNS_DATA\Q.csv');

iu=2;
is=1;
tmax=7200;
tmin=0;
heightlayer1=0.15424;
%% PIEZOMETER P3 PLOT (WF1DUNSAT2)
iu=2;
figure1 = figure('Renderer', 'painters', 'Position', [0 820 800 300],'Color',[1 1 1]);
axes2 = axes('Parent',figure1);
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_hnew(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'--k');
% plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_hsat(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'--k');
hold on
% plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.h_hsat_mean(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'k');
plot(heynsp3.Time,heynsp3.Pressure,'k','LineWidth',2);
hold off
title('Piezometric pressure $\psi_{(s1,u2)}$(m) - Piezometer P3 (Henyns, 2010)','interpreter','latex')
ylabel('Piezometric pressure(m)') 
xlabel('t(s)')
legend({'$\psi_{(s1,u2)}$','Heyns,2010'},'Location','northeast','interpreter','latex')
set(axes2,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')

% %% WATERTABLE (layer 1)
% figure2 = figure('Renderer', 'painters', 'Position', [0 450 800 300],'Color',[1 1 1]);
% axes2 = axes('Parent',figure2);
% hold on
% %ranges of t
% times = [0,10,20,30,40,50,100,200,400,600,800,930];
% for t=times(times>=tmin&times<=tmax)
%     plot(satoutputelements.x0(satoutputelements.is==is&satoutputelements.t==t),satoutputelements.h(satoutputelements.is==is&satoutputelements.t==t),'k');
% end
% times = [950,1000,1500,2000,4000,6000,7200];
% for t=times(times>=tmin&times<=tmax)
%     plot(satoutputelements.x0(satoutputelements.is==is&satoutputelements.t==t),satoutputelements.h(satoutputelements.is==is&satoutputelements.t==t),'-k');
% end
% hold off
% title('Watertables from t='+string(tmin)+'s to t='+string(tmax)+'s. (layer 1)')
% xlabel('Width(m)') 
% ylabel('height(m)') 
% set(axes2,'XMinorGrid','on','YMinorGrid','on');
% xtickformat('%.2f')
%% WF1DSAT: WATERTABLES
iu=1;
is1=1;
is2=2;
tmax=7200;
tmin=0;
figure3 = figure('Renderer', 'painters', 'Position', [0 50 800 300],'Color',[1 1 1]);
axes3 = axes('Parent',figure3);
hold on
%ranges of t
times = [0,10,20,30,40,50,100,200,400,600,800,930];
for t=times(times>=tmin&times<=tmax)
    ly11=plot(satoutputelements.x0(satoutputelements.is==is1&satoutputelements.t==t),satoutputelements.h(satoutputelements.is==is1&satoutputelements.t==t),'k');
end
times = [950,1000,1500,2000,4000,6000,7200];
for t=times(times>=tmin&times<=tmax)
    ly12=plot(satoutputelements.x0(satoutputelements.is==is1&satoutputelements.t==t),satoutputelements.h(satoutputelements.is==is1&satoutputelements.t==t),'--k');
end
times = [0,10,20,30,40,50,100,200,400,600,800,930];
for t=times(times>=tmin&times<=tmax)
    ly21=plot(satoutputelements.x0(satoutputelements.is==is2&satoutputelements.t==t),satoutputelements.h(satoutputelements.is==is2&satoutputelements.t==t)+0.1524,'b');
end
times = [950,1000,1500,2000,4000,6000,7200];
for t=times(times>=tmin&times<=tmax)
    ly22=plot(satoutputelements.x0(satoutputelements.is==is2&satoutputelements.t==t),satoutputelements.h(satoutputelements.is==is2&satoutputelements.t==t)+0.1524,'--b');
end
hold off
title('Watertables from t='+string(tmin)+'s to t='+string(tmax)+'s. (layers 1 and 2, vert module: '+string(iu)+')')
xlabel('Width(m)') 
ylabel('height(m)') 
legend([ly11,ly12,ly21,ly22],{'subballast: t=0-930s','subballast: t=930-7200s','ballast: t=0-930s','ballast: t=930-7200s'},'Location','northeast')
set(axes3,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')

%% WF1DSAT: WATERFLOWS
iu=nan;
is1=1;
is2=2;
tmax=7200;
tmin=0;
figure3 = figure('Renderer', 'painters', 'Position', [0 450 800 300],'Color',[1 1 1]);
axes3 = axes('Parent',figure3);
hold on
%ranges of t
times = [0,10,20,30,40,50,100,200,400,600,800,930];
for t=times(times>=tmin&times<=tmax)
    plot(satoutputelements.x0(satoutputelements.is==is2&satoutputelements.t==t),satoutputelements.qent(satoutputelements.is==is1&satoutputelements.t==t),'k');
end
for t=times(times>=tmin&times<=tmax)
    plot(satoutputelements.x0(satoutputelements.is==is2&satoutputelements.t==t),satoutputelements.dqhordx(satoutputelements.is==is2&satoutputelements.t==t),'g');
end
for t=times(times>=tmin&times<=tmax)
    plot(satoutputelements.x0(satoutputelements.is==is2&satoutputelements.t==t),satoutputelements.dqhordx(satoutputelements.is==is2&satoutputelements.t==t),'--g');
end

times = [950,1000,1500,2000,4000,6000,7200];
for t=times(times>=tmin&times<=tmax)
    plot(satoutputelements.x0(satoutputelements.is==is2&satoutputelements.t==t),satoutputelements.qent(satoutputelements.is==is1&satoutputelements.t==t),'--k');
end
for t=times(times>=tmin&times<=tmax)
    plot(satoutputelements.x0(satoutputelements.is==is2&satoutputelements.t==t),satoutputelements.incvoldt(satoutputelements.is==is2&satoutputelements.t==t),'g');
end
times = [950,1000,1500,2000,4000,6000,7200];
for t=times(times>=tmin&times<=tmax)
    plot(satoutputelements.x0(satoutputelements.is==is2&satoutputelements.t==t),satoutputelements.incvoldt(satoutputelements.is==is2&satoutputelements.t==t),'--g');
end
hold off
title('Waterflows from t='+string(tmin)+'s to t='+string(tmax)+'s. (layers 1 and 2)')
xlabel('Width(m)') 
ylabel('height(m)') 
%legend({'layer 1: t=0-930s','layer 1: t=930-7200s','layer 2: t=0-930s','layer 2: t=930-7200s'},'Location','northeast')
set(axes3,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')



%% WF1DUNSAT(2): qvt vs qvb vs qs  vs t
iu=1;
is=2;
tmax=7200;
tmin=0;
figure1 = figure('Renderer', 'painters', 'Position', [800 50 800 300],'Color',[1 1 1]);
axes2 = axes('Parent',figure1);
% plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_hsat(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'--k');
hold on
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_qvt(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'k');
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_qvb(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'r');
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_Qs_layer(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'b');
hold off
title('WF1DUNSAT: $$q_{vt,u}$$ vs $$q_{vb,u}$$ vs $$Q_{s,u}$$ - iu='+string(iu)+' is= '+string(is),'interpreter','latex')
ylabel('Waterflow(m3/s)') 
xlabel('t(s)')
legend({'$$q_{vt,u}$$','$$q_{vb,u}$$','$$Q_{s,u}$$'},'Location','northeast','interpreter','latex')
set(axes2,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')
% WF1DUNSAT(1): qvt vs qvb vs qs
is=1;
tmax=7200;
tmin=0;
figure1 = figure('Renderer', 'painters', 'Position', [800 50 800 300],'Color',[1 1 1]);
axes2 = axes('Parent',figure1);
% plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_hsat(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'--k');
hold on
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_qvt(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'k');
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_qvb(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'r');
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_Qs_layer(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'b');
hold off
title('WF1DUNSAT: $$q_{vt,u}$$ vs $$q_{vb,u}$$ vs $$Q_{s,u}$$ - iu='+string(iu)+' is= '+string(is),'interpreter','latex')
ylabel('Waterflow(m3/s)') 
xlabel('t(s)')
legend({'$$q_{vt,u}$$','$$q_{vb,u}$$','$$Q_{s,u}$$'},'Location','northeast','interpreter','latex')
set(axes2,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')
%% WF1DUNSAT: incvoldt vs qvt-qvb+qs vs t
iu=2;
is=2;
tmax=7200;
tmin=0;
figure1 = figure('Renderer', 'painters', 'Position', [800 450 800 300],'Color',[1 1 1]);
axes2 = axes('Parent',figure1);
% plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_hsat(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'--k');
hold on
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_incvoldt(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'k');
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_qvt(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax)-modelconstraints.v_qvb(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax)+modelconstraints.v_Qs_layer(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'r');
%plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_Qs_layer(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'b');


%plot(heynsp3.Time,heynsp3.Pressure,'k','LineWidth',2);
hold off
title('WF1DUNSAT: $${\Delta \Theta / dt}_{all}$$ vs $$q_{vt}-q_{vb}-Q_{s,all}$$ - iu='+string(iu)+' is= '+string(is),'interpreter','latex')
ylabel('Waterflow(m3/s)') 
xlabel('t(s)')
legend({'$${\Delta \Theta / dt}_{all}$$','$$\Delta \Theta / dt}_{all}$$'},'Location','northeast','interpreter','latex')
set(axes2,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')

%% hsat_s vs hsat_u (1) vs t
iu=2;
is=1;
tmax=7200;
tmin=0;
figure1 = figure('Renderer', 'painters', 'Position', [800 820 800 300],'Color',[1 1 1]);
axes2 = axes('Parent',figure1);
% plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_hsat(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'--k');
hold on
% plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.h_hsat_mean(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'k');
%h_hsat_mean
is=1;
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.h_hsat_mean(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'k');
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_hsat(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'--k');
is=2;
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.h_hsat_mean(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'b');
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_hsat(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'--b');

%plot(heynsp3.Time,heynsp3.Pressure,'k','LineWidth',2);
hold off
title('Watertable height at P3 (WD1DSAT and WD1DUNSAT)','interpreter','latex')
ylabel('Watertable height(m)') 
xlabel('t(s)')
legend({'$h_{sat,s}$(subballast,WF1DSAT)','$h_{sat,u}$(subballast,WF1DUNSAT)','$h_{sat,s}$(ballast,WF1DSAT)','$h_{sat,u}$(ballast,WF1DUNSAT)'},'Location','northeast','interpreter','latex')
set(axes2,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')

%% WF1DUNSAT: Watercontent vs height
iu=1;
figure2 = figure('Renderer', 'painters', 'Position', [1600 820 400 400],'Color',[1 1 1]);
axes2 = axes('Parent',figure2);
hold on
%ranges of t
times = [0,10,20,30,40,50,100,200,400,600,800,930];
for t=times(times>=tmin&times<=tmax)
    plot(unsatoutputelements.thnew(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'k');
end
times = [950,1000,1500,2000,4000,6000,7200];
for t=times(times>=tmin&times<=tmax)
    plot(unsatoutputelements.thnew(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'-k');
end
hold off
title('WF1DUNSAT('+string(iu)+')Watercontent')
xlabel('Water content(m3/m3)') 
ylabel('height(m)') 
set(axes2,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')

%% WF1DUNSAT: waterflow vs height
iu=3;
figure2 = figure('Renderer', 'painters', 'Position', [1600 420 400 400],'Color',[1 1 1]);
axes2 = axes('Parent',figure2);
hold on

%ranges of t
times = [0,10,20,30,40,50,100,200,400,600,800,930];
for t=times(times>=tmin&times<=tmax)
    plot(unsatoutputelements.qmed(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'k');
%     plot(unsatoutputelements.qs(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'r');
end
times = [950,1000,1500,2000,4000,6000,7200];
for t=times(times>=tmin&times<=tmax)
    plot(unsatoutputelements.qmed(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'--k');
%     plot(unsatoutputelements.qs(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'--r');
end
hold off
title('Waterflow WF1DUNSAT('+string(iu)+')')
xlabel('Waterflow(m3/s)') 
ylabel('height(m)') 
set(axes2,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')
% Waterflow in WF1DUNSAT(2)
iu=2;
figure2 = figure('Renderer', 'painters', 'Position', [1600 420 400 400],'Color',[1 1 1]);
axes2 = axes('Parent',figure2);
hold on
%ranges of t
times = [0,10,20,30,40,50,100,200,400,600,800,930];
for t=times(times>=tmin&times<=tmax)
    plot(unsatoutputelements.qmed(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'k');
%     plot(unsatoutputelements.qs(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'r');
end
times = [950,1000,1500,2000,4000,6000,7200];
for t=times(times>=tmin&times<=tmax)
    plot(unsatoutputelements.qmed(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'--k');
%     plot(unsatoutputelements.qs(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'--r');
end
hold off
title('Waterflow WF1DUNSAT('+string(iu)+')')
xlabel('Waterflow(m3/s)') 
ylabel('height(m)') 
set(axes2,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')

% Waterflow in WF1DUNSAT(1)
iu=1;
figure2 = figure('Renderer', 'painters', 'Position', [1600 420 400 400],'Color',[1 1 1]);
axes2 = axes('Parent',figure2);
hold on
%ranges of t
times = [0,10,20,30,40,50,100,200,400,600,800,930];
for t=times(times>=tmin&times<=tmax)
    plot(unsatoutputelements.qmed(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'k');
%     plot(unsatoutputelements.qs(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'r');
end
times = [950,1000,1500,2000,4000,6000,7200];
for t=times(times>=tmin&times<=tmax)
    plot(unsatoutputelements.qmed(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'--k');
%     plot(unsatoutputelements.qs(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'--r');
end
hold off
title('Waterflow WF1DUNSAT('+string(iu)+')')
xlabel('Waterflow(m3/s)') 
ylabel('height(m)') 
set(axes2,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')



%% WF1DUNSAT: Pore Pressure(m) vs height
iu=3;
figure2 = figure('Renderer', 'painters', 'Position', [1600 20 400 400],'Color',[1 1 1]);
axes2 = axes('Parent',figure2);
hold on
%ranges of t
times = [0,10,20,30,40,50,100,200,400,600,800,930];
for t=times(times>=tmin&times<=tmax)
    plot(unsatoutputelements.hnew(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'k');
end
times = [950,1000,1500,2000,4000,6000,7200];
for t=times(times>=tmin&times<=tmax)
    plot(unsatoutputelements.hnew(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'-k');
end
hold off
title('Pore pressure WF1DUNSAT('+string(iu)+')')
xlabel('Waterflow(m3/s)') 
ylabel('height(m)') 
set(axes2,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')

% Pore Pressure (m) in WF1DUNSAT(2)
iu=2;
figure2 = figure('Renderer', 'painters', 'Position', [1600 20 400 400],'Color',[1 1 1]);
axes2 = axes('Parent',figure2);
hold on
%ranges of t
times = [0,10,20,30,40,50,100,200,400,600,800,930];
for t=times(times>=tmin&times<=tmax)
    plot(unsatoutputelements.hnew(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'k');
end
times = [950,1000,1500,2000,4000,6000,7200];
for t=times(times>=tmin&times<=tmax)
    plot(unsatoutputelements.hnew(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'-k');
end
hold off
title('Pore pressure WF1DUNSAT('+string(iu)+')')
xlabel('Waterflow(m3/s)') 
ylabel('height(m)') 
set(axes2,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')

% Pore Pressure (m) in WF1DUNSAT(1)
iu=1;
figure2 = figure('Renderer', 'painters', 'Position', [1600 20 400 400],'Color',[1 1 1]);
axes2 = axes('Parent',figure2);
hold on
%ranges of t
times = [0,10,20,30,40,50,100,200,400,600,800,930];
for t=times(times>=tmin&times<=tmax)
    plot(unsatoutputelements.hnew(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'k');
end
times = [950,1000,1500,2000,4000,6000,7200];
for t=times(times>=tmin&times<=tmax)
    plot(unsatoutputelements.hnew(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x1(unsatoutputelements.iu==iu&unsatoutputelements.t==t),'-k');
end
hold off
title('Pore pressure WF1DUNSAT('+string(iu)+')')
xlabel('Waterflow(m3/s)') 
ylabel('height(m)') 
set(axes2,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')

%% nrel (1)
iu=2;
is=1;
tmax=7200;
tmin=0;
figure1 = figure('Renderer', 'painters', 'Position', [800 820 800 300],'Color',[1 1 1]);
axes2 = axes('Parent',figure1);
% plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_hsat(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'--k');
hold on
%plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.h_hsat_mean(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'k');
%h_hsat_mean
is=1;
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.nrel(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'k');
is=2;
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.nrel(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'b');

%plot(heynsp3.Time,heynsp3.Pressure,'k','LineWidth',2);
hold off
title('Non-wetted porosity ($$n_{nw}$$)','interpreter','latex')
ylabel('Watertable height(m)') 
xlabel('t(s)')
legend({'nrel(is=1)','nrel(is=2)(m)'},'Location','northeast')
set(axes2,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')

%% WF1DSAT: Last watertable
iu=1;
is1=1;
is2=2;
tmax=7200;
tmin=0;
figure3 = figure('Renderer', 'painters', 'Position', [0 50 800 300],'Color',[1 1 1]);
axes3 = axes('Parent',figure3);
hold on
%ranges of t
times = max(satoutputelements.t);
for t=times(times>=tmin&times<=tmax)
    plot(satoutputelements.x0(satoutputelements.is==is1&satoutputelements.t==t),satoutputelements.h(satoutputelements.is==is1&satoutputelements.t==t),'k');
end
% times = [950,1000,1500,2000,4000,6000,7200];
% for t=times(times>=tmin&times<=tmax)
%     plot(satoutputelements.x0(satoutputelements.is==is1&satoutputelements.t==t),satoutputelements.h(satoutputelements.is==is1&satoutputelements.t==t),'--k');
% end
% times = [0,10,20,30,40,50,100,200,400,600,800,930];
% for t=times(times>=tmin&times<=tmax)
%     plot(satoutputelements.x0(satoutputelements.is==is2&satoutputelements.t==t),satoutputelements.h(satoutputelements.is==is2&satoutputelements.t==t)+0.1524,'b');
% end
% times = [950,1000,1500,2000,4000,6000,7200];
for t=times(times>=tmin&times<=tmax)
    plot(satoutputelements.x0(satoutputelements.is==is2&satoutputelements.t==t),satoutputelements.h(satoutputelements.is==is2&satoutputelements.t==t),'--b');
end
hold off
title('Watertables from t='+string(tmin)+'s to t='+string(tmax)+'s. (layers 1 and 2, vert module: '+string(iu)+')')
xlabel('Width(m)') 
ylabel('height(m)') 
%legend({'layer 1: t=0-930s','layer 1: t=930-7200s','layer 2: t=0-930s','layer 2: t=930-7200s'},'Location','northeast')
set(axes3,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')

%% WF1DUNSAT(2): Waterflows at the boundary (subballast)
iu=1;
is=2;
tmax=7200;
tmin=0;
figure1 = figure('Renderer', 'painters', 'Position', [800 50 800 300],'Color',[1 1 1]);
axes2 = axes('Parent',figure1);
% plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_hsat(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'--k');
hold on
rigthelement=max(satoutputelements.ie);
% qboundballast=satoutputelements.q(satoutputelements.is==2&satoutputelements.ie==rigthelement&satoutputelements.t<=tmax);
% qboundballast_acum=cumsum(qboundballast);
qboundsubballast=satoutputelements.q(satoutputelements.is==2&satoutputelements.ie==rigthelement&satoutputelements.t<=tmax);
qboundsubballast_acum=cumsum(qboundballast);
plot(satoutputelements.t(satoutputelements.is==1&satoutputelements.ie==rigthelement&satoutputelements.t<=tmax),satoutputelements.q(satoutputelements.is==1&satoutputelements.ie==rigthelement&satoutputelements.t<=tmax),'k');
% plot(satoutputelements.t(satoutputelements.is==2&satoutputelements.ie==rigthelement&satoutputelements.t<=tmax),satoutputelements.q(satoutputelements.is==2&satoutputelements.ie==rigthelement&satoutputelements.t<=tmax),'--k');
plot(HeynsQ.t_s(HeynsQ.t_s<=tmax),HeynsQ.Q_subballast(HeynsQ.t_s<=tmax),'b');
plot(HeynsQ.t_s(HeynsQ.t_s<=tmax),HeynsQ.Q_rain(HeynsQ.t_s<=tmax),'r');
hold off
title('Discharge at the seepage boundary (subballast)','interpreter','latex')
ylabel('Waterflow(m3/s)') 
xlabel('t(s)')
legend({'$$Q_{subballast}$$(WF15DSATUNSAT)','$$Q_{subballast}$$ (measured)','$$Q_{rain}$$'},'Location','northeast','interpreter','latex')
set(axes2,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')
%% WF1DUNSAT(2): Waterflows at the boundary (ballast)
iu=1;
is=2;
tmax=7200;
tmin=0;
figure1 = figure('Renderer', 'painters', 'Position', [800 50 800 300],'Color',[1 1 1]);
axes2 = axes('Parent',figure1);
% plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_hsat(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),'--k');
hold on
rigthelement=max(satoutputelements.ie);
qboundballast=satoutputelements.q(satoutputelements.is==2&satoutputelements.ie==rigthelement&satoutputelements.t<=tmax);
qboundballast_acum=cumsum(qboundballast);
% qboundsubballast=satoutputelements.q(satoutputelements.is==2&satoutputelements.ie==rigthelement&satoutputelements.t<=tmax);
% qboundsubballast_acum=cumsum(qboundballast);
% plot(satoutputelements.t(satoutputelements.is==1&satoutputelements.ie==rigthelement&satoutputelements.t<=tmax),satoutputelements.q(satoutputelements.is==1&satoutputelements.ie==rigthelement&satoutputelements.t<=tmax),'k');
plot(satoutputelements.t(satoutputelements.is==2&satoutputelements.ie==rigthelement&satoutputelements.t<=tmax),satoutputelements.q(satoutputelements.is==2&satoutputelements.ie==rigthelement&satoutputelements.t<=tmax),'--k');
plot(HeynsQ.t_s(HeynsQ.t_s<=tmax),HeynsQ.Q_ballast(HeynsQ.t_s<=tmax),'b');
% plot(HeynsQ.t_s(HeynsQ.t_s<=tmax),HeynsQ.Q_subballast(HeynsQ.t_s<=tmax),'b');
plot(HeynsQ.t_s(HeynsQ.t_s<=tmax),HeynsQ.Q_rain(HeynsQ.t_s<=tmax),'r');
hold off
title('Discharge at the seepage boundary (ballast)','interpreter','latex')
ylabel('Waterflow(m3/s)') 
xlabel('t(s)')
legend({'$$Q_{ballast}$$(WF15DSATUNSAT)','$$Q_{ballast}$$ (measured)','$$Q_{rain}$$'},'Location','northeast','interpreter','latex')
set(axes2,'XMinorGrid','on','YMinorGrid','on');
xtickformat('%.2f')