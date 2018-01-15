%% SAMPLE FORCE INPUT #1: FREE VIBRATION
% p_v=zeros(nn,1,(length(t))); %a zero vector matching the size of t

%% SAMPLE FORCE INPUT #2:5-Story Chimney with step force applied at top
% po=[0;0;0;0;0]; %initial force -kips [assume 0 at t=0 of step force]
% p=[0;0;0;0;1000]; %kips; force at every node--@t>0 of step force
% %^as it's written, the top story will have a step force of 1000 kips
% %applied to it.  
% 
% p_v=zeros(nn,1,(length(t))); %create a stack of column vectors
% p_v(:,1,1)=po; %assign its initial column vectors the initial force vector
% for i=2:length(t)
%     p_v(:,1,i)=p; %assign a steady force of p vector for the subsequent 
%     %force column vectors 
% end

%% SAMPLE FORCE INPUT #3: El_Centro Ground Motion
%[M]{a}+[C]{v}+[K]{u}=-[M]{i}{ag}
run('El_Centro_Ground_Motion')
b=a'; 
zz=b(:); %c turns data into neat column vector;
zz(length(zz))=[]; %last value in motion data 0.000 was added to complete the square matrix
acc=[0;zz]; %@ t=0, acceleration is 0, motion data starts at t=0.02;
el_t=0:0.02:31.18; %the time vector goes from 0 to 31.18 with time step of .02s;

x=el_t'; %El Centro Time vector converts to column vector
y=acc; %El Centro acc vector already column vector.
xti=t'; %convert user input time vector into column vector
yi=interp1q(x,y,xti); %uses linear interpolation to get data points spaced dt apart.  
yi(isnan(yi))=0; %values after the last data in ground motion acceleration=0.  
g=32.2*12; %gravity in inches/s^2
ag=yi*g; %ground motion acc. a(t);

influ=ones(nn,1); %since only a chimney is being modeled.  
p_v=zeros(nn,1,(length(t))); %create a stack of column vectors
for i=1:length(t)
    p_v(:,1,i)=-m*influ*ag(i);
end

%% SAMPLE FORCE INPUT #4: #-story chimney w.harmonic function at top DOF
% ptr=tr; %forcing function time range [s]
%     %set as equal to time range for now...
% p_amp=1000; %forcing amplitude [kips]
% wf=pi/0.6; %forced frequency [rad/s]
% p_v=zeros(nn,1,(length(t))); %create a stack of column vectors
% for i=1:length(t)
%     if t(i)<=ptr %during the interval the external force is applied...
%     p_v(nn,1,i)=p_amp*sin(wf*t(i)); 
%     %^^apply to the uppermost dof (nn) this harmonic forcing function.  
%     else
%         p_v(i)=0; %outside of the applied force range, the force remains zero.
%     end
% end


 