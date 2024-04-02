% constructSimInput writes the variables needed for the simulation to run
% to its workspace. author: ishan
function constructSimInput(MOI, l, m, g, inputLimits, throttleConsts, betaInputDelay, gammaInputDelay, discretize_interval, ...
   x_set, u_set, t_set, Kset, tSegs, startTime, stopTime, q1, x1)
   mdlWks = get_param('Model_Mk2','ModelWorkspace'); % model workspace
   clear(mdlWks); % start from scratch

   % todo: add validation
   assignin(mdlWks,'MOI',MOI);
   assignin(mdlWks,'l',l);
   assignin(mdlWks,'m',m);
   assignin(mdlWks,'g',g);
   assignin(mdlWks,'inputLimits',inputLimits);
   assignin(mdlWks,'throttleConsts',throttleConsts);
   assignin(mdlWks,'betaInputDelay',betaInputDelay);
   assignin(mdlWks,'gammaInputDelay',gammaInputDelay);
   assignin(mdlWks,'discretize_interval',discretize_interval);
   assignin(mdlWks,'x_set',x_set);
   assignin(mdlWks,'u_set',u_set);
   assignin(mdlWks,'t_set',t_set);
   assignin(mdlWks,'Kset',Kset);
   assignin(mdlWks,'tSegs',tSegs);
   assignin(mdlWks,'startTime',startTime);
   assignin(mdlWks,'stopTime',stopTime);

   % these ones probably shouldn't be needed tbh. we should change this
   assignin(mdlWks,'q1',q1);
   assignin(mdlWks,'x1',x1);

   % saveToSource(mdlWks)

   save_system('Model_Mk2')

end