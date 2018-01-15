%% Runge-Kutta program 
run('ForcingFunction')


%% initialization matrices
du1=zeros(1,length(t));dv1=zeros(1,length(t));du2=zeros(1,length(t));
dv2=zeros(1,length(t));du3=zeros(1,length(t));dv3=zeros(1,length(t));
du4=zeros(1,length(t));dv4=zeros(1,length(t));du=zeros(1,length(t));
dv=zeros(1,length(t));u=zeros(1,length(t));v=zeros(1,length(t));
%% Runge-Kutta 4th Order ODE Iterations
%initial conditions
v(1)=v0; %set 1st pt. of velocity vector as initial velocity
u(1)=u0; %set 1st pt. of displacement vector as initial displacement

rkc=c/m; %normalize damping coefficient by mass.
rkk=k/m; %normalize stiffness coefficient by m

for i=2:length(t) %for every time step dt...
du1(i)=dt*v(i-1); %Eqn 1.7.2.1
dv1(i)=dt*((p(i-1))/m-rkc*v(i-1)-rkk*u(i-1)); %Eqn 1.7.2.2
        
du2(i)=dt*(v(i-1)+(dv1(i)/2)); %Eqn 1.7.2.3
dv2(i)=dt*((.5*(p(i)+p(i-1)))/m-rkc*(v(i-1)+dv1(i)/2)-rkk*(u(i-1)+du1(i)/2));
%^Eqn 1.7.2.4
  
du3(i)=dt*(v(i-1)+dv2(i)/2); %Eqn 1.7.2.5
dv3(i)=dt*((.5*(p(i)+p(i-1)))/m-rkc*(v(i-1)+dv2(i)/2)-rkk*(u(i-1)+du2(i)/2));
%^Eqn 1.7.2.6

du4(i)=dt*(v(i-1)+dv3(i)); %Eqn 1.7.2.7
dv4(i)=dt*(((p(i)))/m-rkc*(v(i-1)+dv3(i))-rkk*(u(i-1)+du3(i))); %Eqn 1.7.2.8
 
du(i)=(du1(i)+2*du2(i)+2*du3(i)+du4(i))/6; %Eqn 1.7.2.9
dv(i)=(dv1(i)+2*dv2(i)+2*dv3(i)+dv4(i))/6; %Eqn 1.7.2.10
       
 u(i)=u(i-1)+du(i); %Eqn 1.7.2.11
 v(i)=v(i-1)+dv(i); %Eqn 1.7.2.12
end 
    %% ^^see pseudocode off YouTube https://www.youtube.com/watch?v=smfX0Jt_f0I
%9:11 into video
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

%% Displacement by Runge Kutta vs. time
subplot(3,1,2) %the middle graph
plot(t,u,'b'); %plot it in blue
xlabel('t, [seconds]','FontSize',13)
ylabel('u,[inches]','FontSize',13)
title('Displacement by Runge-Kutta','FontSize',16)
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
  %h=traces along spring graph
    set(h,'XData',u(i)); %for the 1D plot
    [xs,ys]=spring(min(u)*1.1,0,u(i),0);
      plot(xs,ys,'LineWidth',2,'Color','green') %plot a green spring
    pause(dt); %animate in real time
    if t(i)<t(end)
     delete(findobj('Color','green')) %a way to animate the spring
     %delete everything that's green in the figure (only the spring)
     %for every frame of the loop.  
    set(gcf,'Renderer','OpenGL')
    end
    i=i+1;
end
toc





