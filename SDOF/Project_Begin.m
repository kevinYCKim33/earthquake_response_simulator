%% Run this program after user_inputs all values in User_Input.m and 
%enters forcing function in ForcingFunction.m  
clc;
clear all;
close all;

run('User_Input')

t=(0:dt:tr);%time vector [s]; from 0 to tr by time step dt.

wn=sqrt(k/m); %natural frequency (rad/s)
c=2*m*wn*xi; %damping coefficient (kip*s/in)
Tn= 2*pi/wn; %natural period (s)


switch lower(Analysis_Type) %turns the user-inputs into lowercase letters

case 'cdm' %for Central Difference Method Case...
    if Ground_Motion>1 %if user is using smooth forcing function
        run('Central_Difference_Method') %run this program
    else
        run('Central_Difference_Method_Ground_Motion') %run this for 
        %running El_Centro Ground Motion
    end

case 'newmark' %for Newmark's Method Case
     if Ground_Motion>1 %if user is using smooth forcing function
        run('Newmark_Method') %% run this program 
     else
        run('Newmark_Method_Ground_Motion') %% run this program
        %for running El_Centro Ground Motion
     end
     
case 'rungekutta' %for the Runge-Kutta Case
     if Ground_Motion>1 %if user is using smooth forcing function
        run('Runge_Kutta_Method') %% run this program 
     else
        run('Runge_Kutta_Method_Ground_Motion') %% run this program
        %for running El_Centro Ground Motion
     end
end