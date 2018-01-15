%% User_input.m--attain output by running Project_Begin.m

%% structrual properties
m=0.2533; %kip-sec^2/in %mass of SDOF structure
k=10; %kips/in %stiffness of SDOF structure
xi=0.05; %damping ratio

%% initial conditions
u0=0; %initial displacement (in)
v0=0; %initial velocity (in/s)

%% numerical integration
dt=.1; %time increment (s)
tr=10; %time range (s)

Analysis_Type='CDM'; %Enter 'CDM','Newmark',or 'RungeKutta' 
NM_Type='N/A'; %Enter 'average','linear', or'N/A' for non-Newmark analysis

Ground_Motion=2;  %Enter 1 if feeding recorded ground motion into prog.
%ForcingFunction; %Enter 2 if inputting theoretical forcing functions.

%% Head into ForcingFunction.m to assign a forcing fcn. or ground motion



