set_param('moore_machine','SaveFinalState','on','FinalStateName',...
'myOperPoint','SaveOperatingPoint','on');
simOut = sim('moore_machine','StartTime','0','StopTime','9')
myOperPoint = simOut.myOperPoint

set_param('moore_machine','LoadInitialState','on','InitialState',...
'myOperPoint');
myOperPoint = simOut.myOperPoint
simOut = sim('moore_machine','StartTime','9','StopTime','18')

set_param('moore_machine','LoadInitialState','off','InitialState',...
'myOperPoint');
myOperPoint = simOut.myOperPoint
simOut = sim('moore_machine','StartTime','9','StopTime','18')
