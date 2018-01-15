 %%Newmark's Method
run ('ForcingFunction'); %Forcing Function creates a p vector spaced dt apart

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

u(1)=u0; %replace the first u entry as initial displacement
v(1)=v0;
acc_n(1)=a0; 

%% Table 5.4.2 (213) 2.0 Calculations for each time step, i
for i=1:length(t)-1
    dp(i)=p(i+1)-p(i);

    dph(i)=dp(i)+a*v(i)+b*acc_n(i);
    du(i)=dph(i)/kh;
    dv(i)=(g/B/dt)*du(i)-(g/B)*v(i)+dt*(1-(g/2/B))*acc_n(i);
    da(i)=1/(B*dt^2)*du(i)-(1/(B*dt))*v(i)-(1/(2*B))*acc_n(i);
    
    u(i+1)=u(i)+du(i);
    v(i+1)=v(i)+dv(i);
    acc_n(i+1)=acc_n(i)+da(i);
end

%% plot the displacement results
figure('units','normalized','outerposition',[0 0 1 1]) %full screens the figures

%% El_Centro Acceleration Time History
subplot(2,2,1) 
plot(el_t,acc,'r'); %plot acceleration history  
hold on
gg=plot(xi,yi,'y.'); %plot the interpolated values of ground motion with yellow dots.
gh=plot(xi(1),yi(1),'bx','MarkerSize',20,'MarkerFaceColor','b','LineWidth',2); 
    %^^put a crosshair later to be animated at initial conditions
    grid
    xlabel('time (s)','FontSize',13)
    ylabel('Acceleration (g)','FontSize',13)
    title('El Centro Acceleration Time History','FontSize',16)
    xlim([0,tr]);

%% Displacement by Newmark vs. time
subplot(2,2,3) 
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

p=plot(t(1),u(1),'rx','MarkerSize',20,'MarkerFaceColor','b','LineWidth',2);
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



