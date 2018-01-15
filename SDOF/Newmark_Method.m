 %%Newmark's Method
run ('ForcingFunction'); %needs the forcing function p from it.  

%% takes the user inputs of NM by average or linear acc. method, to
%determine the gamma and beta values
switch lower(NM_Type)
       case 'average' 
           B=.25;
           g=.5;
       case 'linear'
           B=1/6;
           g=.5;
            if dt/Tn>.551
                error('time step too large')
            end
end
%% initial calculations Newmark. look at page 213, Table 5.4.2
a0=(p(1)-c*v0-k*u0)/m; %initial acceleration 1.1
kh=k+(g/(B*dt))*c+(m/(B*dt^2)); %Step 1.3
a=(m/(B*dt))+(g*c/B); b=m/(2*B)+dt*(g/(2*B)-1)*c; %Step 1.4

%% Calculations for each time step i
u=zeros(1,length(t)); %initialize a displacement vector same size as t vec.
v=zeros(1,length(t)); % initialize velocity vector
acc_n=zeros(1,length(t)); %initialize accelearation vector

%initializing various intermediate variables
dp=zeros(1,length(t));dph=zeros(1,length(t));du=zeros(1,length(t)); 
dv=zeros(1,length(t));da=zeros(1,length(t));  

%% Initial conditions
u(1)=u0; %replace the first u entry as initial displacement
v(1)=v0; %replace the first velocity entry as initial velocity
acc_n(1)=a0; %replace the first acceleration entry as initial acceleration
%^^where a0 is Eqn 1.6.2.1

% Iterative Steps
for i=1:length(t)-1
    dp(i)=p(i+1)-p(i); %Eqn 1.6.2.5
    dph(i)=dp(i)+a*v(i)+b*acc_n(i);%Eqn 1.6.2.6
    du(i)=dph(i)/kh; %Eqn 1.6.2.7
    dv(i)=(g/B/dt)*du(i)-(g/B)*v(i)+dt*(1-(g/2/B))*acc_n(i);%1.6.2.8
    da(i)=1/(B*dt^2)*du(i)-(1/(B*dt))*v(i)-(1/(2*B))*acc_n(i);%1.6.2.9
    
    u(i+1)=u(i)+du(i);%1.6.2.10
    v(i+1)=v(i)+dv(i);%1.6.2.11
    acc_n(i+1)=acc_n(i)+da(i);%1.6.2.12
end

%% plot the displacement results
figure('units','normalized','outerposition',[0 0 1 1]) %full screens the figures

%% Forcing Function plot
subplot(3,1,1);%the top graph
plot(t,p,'r'); %plot forcing function in red.
xlabel('t, [seconds]','FontSize',13)
ylabel('p,[kips]','FontSize',13)
title('External Force vs. Time','FontSize',16)
grid on
hold on
set(gcf,'Renderer','OpenGL'); %use computer graphic card if any,

f=plot(t(1),p(1),'bx','MarkerSize',10,'MarkerFaceColor','b');
% Put blue crosshair at initial time later to be animated.

%% Displacement by Newmark vs. time
subplot(3,1,2) 
plot(t,u,'b'); %plot it in blue
xlabel('t, [seconds]','FontSize',13)
ylabel('u,[inches]','FontSize',13)
if B>1/6
    title('Displacement by Newmark-Avg. Acc','FontSize',15)
else
    title('Displacement by Newmark-Linear Acc','FontSize',15)  
end
grid on
hold on
set(gcf,'Renderer','OpenGL'); %use computer graphic card if any,

pp=plot(t(1),u(1),'rx','MarkerSize',10,'MarkerFaceColor','r');
% again put crosshair at initial conditions later to be animated.

%% Spring Simulation  
subplot(3,1,3)
h=plot(u(1),0,'rs','MarkerSize',30,'MarkerFaceColor','r');
%plot a big red square at initial position u0 to act as the mass.
hold on
[xs,ys]=spring(min(u)*1.1,0,0,0,round(max(u)*10),1,0.1);
%^^a downloaded function that draws a spring positioned at x0,y0,
%and x1,y1 on other end, and assigns a number of coil.
%see function spring.m for more detail.  
zz=plot(xs,ys,'LineWidth',2,'Color','green'); %plot the spring in green.  

set(gcf,'Renderer','OpenGL'); %use graphic card of computer to animate.  

xlim([min(u)*1.1,max(u)*1.1]); %set x-window range to be 10% larger than 
%max and min displacement the mass will undergo...
ylim([-1.0,1.0]); %arbitrary y range
set(gca,'YTick',[]); %omit y-axis label, only concerned with 1D animation.  
grid on
xlabel('u,displacement (in)','FontSize',13);
title('1D Spring Motion','FontSize',16);


%% Animation Loop
i=1;
tic
while i<=length(u)
  %f=traces along forcing graph
    set(f,'XData',t(i));
    set(f,'YData',p(i));
  %pp=traces along displacement graph
    set(pp,'XData',t(i));
    set(pp,'YData',u(i));
  %h=traces along structural simulation graph
    set(h,'XData',u(i)); %for the 1D plot
    [xs,ys]=spring(min(u)*1.1,0,u(i),0);
      plot(xs,ys,'LineWidth',2,'Color','green') %plot a green spring
    pause(dt) %animate in real time 
    if t(i)<t(end)
     delete(findobj('Color','green')) %a way to animate the spring
     %delete everything that's green in the figure (only the spring)
     %for every frame of the loop.  
    set(gcf,'Renderer','OpenGL')
    end
    i=i+1;
end
toc


