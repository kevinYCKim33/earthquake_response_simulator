%% Central Difference Method Case (Sec. 5.3 for proof)
run ('ForcingFunction'); %Forcing Function of El Centro Ground Motion
%% initial calculations for CDM
a0=(p(1)-(c*v0)-(k*u0))/m; %u double dot; initial acceleration (in/s^2)
u_1=u0-dt*v0+(dt)^2*a0/2; %eqn 1.1 in Example 5.2 (209)
kh=(m/dt^2)+(c/(2*dt)); %khat eqn 1.3 Example 5.2 
a= (m/dt^2)-(c/(2*dt)); %a, eqn 1.4 Example 5.2 
b= k-((2*m)/(dt)^2);  %b, eqpn 1.5 Example 5.2

%% solution blow up?
if (dt/Tn)>(1/pi) %eqn (5.3.11)
    error('time step too large for CDM to be any accurate') 
end

%%calculations for time step i
u=zeros(1,length(t)); %initialize a displacement vector same size as t vec.
phat=zeros(1,length(t)); %initialize a p vector same size as t vec.
u(1)=u0; %replace the first u entry as user input initial displacement
for i=2:length(t)
    phat(i)=p(i)-a*u(i-1)-b*u(i); %Equation 1.5.1.4
    u(i+1)=phat(i)/kh;            %Equation 1.5.1.3 
end
 u(length(u))=[]; %take out the last entry in u vector since u 
%vector is one size bigger than it needs to be.


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

f=plot(t(1),p(1),'bx','MarkerSize',20,'MarkerFaceColor','b');
% Put blue crosshair at initial time later to be animated.

%% Displacement by CDM vs. time
subplot(3,1,2) %the middle graph
plot(t,u,'b'); %plot it in blue
xlabel('t, [seconds]','FontSize',13)
ylabel('u,[inches]','FontSize',13)
title('Displacement by Central Difference Method','FontSize',16)
grid on
hold on
set(gcf,'Renderer','OpenGL'); %use computer graphic card if any,

pp=plot(t(1),u(1),'rx','MarkerSize',20,'MarkerFaceColor','r');
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
tic %measures code efficiency
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
    pause(dt)
    if t(i)<t(end)
     delete(findobj('Color','green')) %a way to animate the spring
     %delete everything that's green in the figure (only the spring)
     %for every frame of the loop.  
    set(gcf,'Renderer','OpenGL')
    end
    i=i+1;
end
toc 



