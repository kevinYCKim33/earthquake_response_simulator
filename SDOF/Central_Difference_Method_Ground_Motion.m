%% Central Difference Method Case (210 for proof)
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
u(1)=u0; %replace the first u entry as initial displacement

for i=2:length(t)
    phat(i)=p(i)-a*u(i-1)-b*u(i); %equation 2.1 Example 5.2
    u(i+1)=phat(i)/kh; %equation 2.2 Example 5.2
end
 u(length(u))=[]; %take out the last entry in u vector since u 
%vector is one size bigger than it needs to be.
%% plot the displacement results
figure('units','normalized','outerposition',[0 0 1 1]) %full screens the figures

%% El_Centro Acceleration Time History
subplot(2,2,1) 
plot(el_t,acc); %plot acceleration history  
hold on
gg=plot(xi,yi,'y.'); %plot the interpolated values of ground motion with yellow dots.
gh=plot(xi(1),yi(1),'rx','MarkerSize',10,'MarkerFaceColor','r'); 
    %^^put a crosshair later to be animated at initial conditions
    grid
    xlabel('time (s)','FontSize',13)
    ylabel('Acceleration (g)','FontSize',13)
    title('El Centro Acceleration Time History','FontSize',16)

%% Displacement by CDM vs. time
subplot(2,2,3) 
plot(t,u,'b'); %plot it in blue
xlabel('t, [seconds]','FontSize',13)
ylabel('u,[inches]','FontSize',13)
title('Displacement by Central Difference Method','FontSize',16)
grid on
hold on
set(gcf,'Renderer','OpenGL'); %use computer graphic card if any,

p=plot(t(1),u(1),'rx','MarkerSize',10,'MarkerFaceColor','r');
% again put crosshair at initial conditions later to be animated.

%% Structure Simulation  
subplot(2,2,[2 4]) %takes up the right half of the screen
h=max(u)*50; %height of arbitary rod set 50x the max value of horizontal disp.
j=plot(u(1),h,'ro','MarkerSize',20,'MarkerFaceColor','r'); 
% ^plot a big red circle at the top of the arbitrary column
hold on
xlim([min(u)*2.5,max(u)*2.5]); %ranges from minimum value and max val. *1.1
ylim([0,h*1.5]); %y range window
set(gca,'YTick',[]); %makes the y-axis label go away.  
grid on
xlabel('u,displacement (in)','FontSize',13);
title('SDOF Disp w. Mass m, Stiffness k','FontSize',16);

%% Animation Loop
i=1;
tic
while i<=length(u)
  %gh=traces along acceleration graph graph
    set(gh,'XData',xi(i));
    set(gh,'YData',yi(i));
  %p=traces along displacement graph
    set(p,'XData',t(i));
    set(p,'YData',u(i));
  %j=traces along structural simulation graph
    set(j,'XData',u(i)); %for the 1D plot
    s01=line([0,u(i)],[0,h],'Color','k','LineWidth',4); %draw black line from
    %base to the mass.  
    pause(dt); %animate in real time
    if t(i)<t(end)
    delete(findobj(gca,'Color','black')); %Trick to animating the rod.  
    set(gcf,'Renderer','OpenGL')
    end
    i=i+1;
end
toc



