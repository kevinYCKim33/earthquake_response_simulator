%% Runge-Kutta program 
run('ForcingFunction')
rkc=c/m; %normalize damping coefficient by mass.
rkk=k/m; %normalize stiffness coefficient by m

%% initialization matrices
du1=zeros(1,length(t));dv1=zeros(1,length(t));du2=zeros(1,length(t));
dv2=zeros(1,length(t));du3=zeros(1,length(t));dv3=zeros(1,length(t));
du4=zeros(1,length(t));dv4=zeros(1,length(t));du=zeros(1,length(t));
dv=zeros(1,length(t));u=zeros(1,length(t));v=zeros(1,length(t));

%initial conditions
v(1)=v0;
u(1)=u0;

%% see pseudocode off YouTube https://www.youtube.com/watch?v=smfX0Jt_f0I
%9:11 into video
    for i=2:length(t)
        du1(i)=dt*v(i-1);
        dv1(i)=dt*((p(i-1))/m-rkc*v(i-1)-rkk*u(i-1));
        
        du2(i)=dt*(v(i-1)+(dv1(i)/2));
        dv2(i)=dt*((.5*(p(i)+p(i-1)))/m-rkc*(v(i-1)+dv1(i)/2)-rkk*(u(i-1)+du1(i)/2));
        %^^linear interpolation half step
        
        du3(i)=dt*(v(i-1)+dv2(i)/2);
        dv3(i)=dt*((.5*(p(i)+p(i-1)))/m-rkc*(v(i-1)+dv2(i)/2)-rkk*(u(i-1)+du2(i)/2));
        %^^linear interpolation half step
        
        du4(i)=dt*(v(i-1)+dv3(i));
        dv4(i)=dt*(((p(i)))/m-rkc*(v(i-1)+dv3(i))-rkk*(u(i-1)+du3(i)));
        %^^linear interpolation FULL step
        
        du(i)=(du1(i)+2*du2(i)+2*du3(i)+du4(i))/6;
        dv(i)=(dv1(i)+2*dv2(i)+2*dv3(i)+dv4(i))/6;
       
        u(i)=u(i-1)+du(i);
        v(i)=v(i-1)+dv(i);
    end
%% plot the displacement results
figure('units','normalized','outerposition',[0 0 1 1])%fullscreens figures

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

%% Displacement by Runge-Kutta vs. time
subplot(2,2,3) 
plot(t,u,'b'); %plot it in blue
xlabel('t, [seconds]','FontSize',13)
ylabel('u,[inches]','FontSize',13)
title('Displacement by Runge-Kutta','FontSize',16)
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

