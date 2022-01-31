%Plot: iu in element 1 for t<tmax
% unsatoutputelements = import_csv('unsat_outputs_elements.csv');
% satoutputelements   = import_csv('sat_outputs_elements.csv');
% modelconstraints    = import_csv('model_constraints.csv');
inputfilename = 'model_inputs.wfuinp';
[p,f,e]=fileparts(inputfilename);
inputfilename=fullfile(p,f);
unsatoutputelements = readtable(strcat(inputfilename,'.outelmu'),'FileType','delimitedtext');
satoutputelements = readtable(strcat(inputfilename,'.outelms'),'FileType','delimitedtext');
modelconstraints = readtable(strcat(inputfilename,'.outcons'),'FileType','delimitedtext');

% unsatoutputelements = readtable('unsat_outputs_elements.csv');
% satoutputelements   = readtable('sat_outputs_elements.csv');
% modelconstraints    = readtable('model_constraints.csv');
iu=2;
nelemv=20;
nelemh=1;
tmin = 0;
tmax = 7200;
is=1;
times = unique(satoutputelements.t)';
close all
figure
plot(unsatoutputelements.t(unsatoutputelements.nelem==nelemv&unsatoutputelements.iu==iu&unsatoutputelements.t<=tmax),unsatoutputelements.hnew(unsatoutputelements.nelem==nelemv&unsatoutputelements.iu==iu&unsatoutputelements.t<=tmax));
title('Piezometer P1(1.5D model vs Heyns,2000)')
xlabel('t(s)') 
ylabel('pore pressure(mH2O)')
figure
for t=times(times>=tmin&times<=tmax)
plot(unsatoutputelements.thnew(unsatoutputelements.iu==iu&unsatoutputelements.t==t),unsatoutputelements.x0(unsatoutputelements.iu==iu&unsatoutputelements.t==t));
hold on
end
%% WATERFLOW: WD1DSat, ELEMENT: 
title('Water content in left vertical unsaturated element from t='+string(tmin)+'s to t='+string(tmax)+'s.')
xlabel('theta (water content) (m3/m3)') 
ylabel('height (m)')
figure
plot(satoutputelements.t(satoutputelements.ie==max(satoutputelements.ie)&satoutputelements.is==is&satoutputelements.t>=tmin&satoutputelements.t<=tmax),satoutputelements.q_all(satoutputelements.ie==max(satoutputelements.ie)&satoutputelements.is==is&satoutputelements.t>=tmin&satoutputelements.t<=tmax));
hold on
plot(satoutputelements.t(satoutputelements.ie==max(satoutputelements.ie)&satoutputelements.is==is&satoutputelements.t>=tmin&satoutputelements.t<=tmax),satoutputelements.q(satoutputelements.ie==max(satoutputelements.ie)&satoutputelements.is==is&satoutputelements.t>=tmin&satoutputelements.t<=tmax));
title('Waterflows leaving(1.5D model vs Heyns,2000)')
xlabel('t(s)') 
ylabel('q(m3/s)')
legend({'total q','1 layer 1'},'Location','northeast')
hold off
figure
for t=times(times>=tmin&times<=tmax)
    plot(satoutputelements.x0(satoutputelements.is==is&satoutputelements.t==t),satoutputelements.h(satoutputelements.is==is&satoutputelements.t==t));
    hold on
end
hold off
title('Watertables from t='+string(tmin)+'s to t='+string(tmax)+'s.')
xlabel('Width(m)') 
ylabel('height(m)') 

figure
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_hsat(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax));
hold on
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.h_hsat_mean(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax));
hold off
title('hsat in unsat and in sat (m)')
ylabel('Height(m)') 
xlabel('t(s)')
legend({'hunsat','hsat'},'Location','northeast')

figure
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_qvtb(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax));
hold on
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_incvoldt(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax));
hold on
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.h_incvoldt_mean(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax));
hold on
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_Qs_layer(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax));
hold on
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.h_dqhordx_mean(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax));
hold on
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.v_Qs_all(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax));
hold on
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.h_dqhordx_all_mean(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax));
hold off
title('Constrain: Sat: '+string(is)+' Unsat: '+string(iu)+' Flows in each timestep')
ylabel('q(m3/s)') 
xlabel('t(s)')
legend({'v:qver','v:incvoldt','h:incvoldt','v:qnewmann','h:dqhordx','v:qnewmann_{all}','h:dqhordx_{all}'},'Location','northeast')

figure
plot(modelconstraints.t(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax),modelconstraints.nrel(modelconstraints.is==is&modelconstraints.iu==iu&modelconstraints.t<=tmax));
title('nrel in time')
ylabel('nrel') 
xlabel('t(s)')

figure
plot(satoutputelements.t(satoutputelements.ie==nelemh&satoutputelements.is==is&satoutputelements.t>=tmin&satoutputelements.t<=tmax),satoutputelements.qent(satoutputelements.ie==nelemh&satoutputelements.is==is&satoutputelements.t>=tmin&satoutputelements.t<=tmax));
hold on
plot(satoutputelements.t(satoutputelements.ie==nelemh&satoutputelements.is==is&satoutputelements.t>=tmin&satoutputelements.t<=tmax),satoutputelements.incvoldt(satoutputelements.ie==nelemh&satoutputelements.is==is&satoutputelements.t>=tmin&satoutputelements.t<=tmax));
hold on
plot(satoutputelements.t(satoutputelements.ie==nelemh&satoutputelements.is==is&satoutputelements.t>=tmin&satoutputelements.t<=tmax),satoutputelements.dqhordx_from_incvoldt_all(satoutputelements.ie==nelemh&satoutputelements.is==is&satoutputelements.t>=tmin&satoutputelements.t<=tmax));
hold on
plot(satoutputelements.t(satoutputelements.ie==nelemh&satoutputelements.is==is&satoutputelements.t>=tmin&satoutputelements.t<=tmax),satoutputelements.dqhordx_from_incvoldt(satoutputelements.ie==nelemh&satoutputelements.is==is&satoutputelements.t>=tmin&satoutputelements.t<=tmax));
hold on
title('Sat model:'+string(is)+' - Waterflows on element '+string(nelemh))
xlabel('t(s)') 
ylabel('q(m3/s)')
legend({'qent','incvoldt','dqhordx from incvol all','dqhordx from incvol'},'Location','northeast')
hold off
% legend({'y = sin(x)','y = cos(x)'},'Location','southwest')


