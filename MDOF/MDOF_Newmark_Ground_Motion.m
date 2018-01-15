run('MDOF_ForcingFunction'); %user defined fun

%% MDOF Newmark acceleration: Section 2.5.1 of Report
%Also mentioned in Chopra Table 15.2.2 pg. 650

%coefficients: assumes 'average'; soln blows up with modes bigger than 2
%when 'linear' option chosen
beta=.25;
y=.5; %gamma

check=size(phi); %gets the x by y dimensions of PHI matrix
w=check(2); %i'm just interested in the width

%initial calculations
%1.1
    qo_v=zeros(1,length(nn)); 
    dqo_v=zeros(1,length(nn));
for i=1:w %for every column of phi vector; (for each mode shape)
    qo_v(i)=(phi(:,i)'*m*u0)/(phi(:,i)'*m*phi(:,i)); 
    %^^ Eqn (2.5.1.1):get initial modal coordinate qo
    dqo_v(i)=(phi(:,i)'*m*v0)/(phi(:,i)'*m*phi(:,i)); 
    %^^ Eqn (2.5.1.1):get initial modal "velocity" dqo
end
    qo_v=qo_v'; %transpose as column vector
    dqo_v=dqo_v'; %transpose as column vector
% %1.2
Po=phi'*p_v(:,1,1); %get initial P matrix normalized by phi.
% ^Eqn (2.5.1.2)

% %1.3: solve for double derivative of qo.
RHS=Po-C*dqo_v-K*qo_v; 
aqo=M^(-1)*RHS;
% ^Eqn (2.5.1.3)

% %1.4 Select dt, and tr from User_Input_MDOF.

% %1.5 Get capital K hat
Kh=K+(y/(beta*dt))*C+(1/(beta*(dt)^2))*M; %Eqn 2.5.1.4
% 
%1.6 get a and b diagonal matrices
a=(1/(beta*dt))*M+(y/beta)*C; %Eqn 2.5.1.5
b=(1/(2*beta))*M+dt*(y/(2*beta)-1)*C;%Eqn 2.5.1.6
% 
%2.0 time step calculations

%% Create Capital P Vector: phi^T*p
P_v=zeros(nn,1,length(t));
for i=1:length(t)
    P_v(:,:,i)=phi'*p_v(:,:,i); %Eqn 2.5.1.7
end

%%initialize q(0), q'(0),q"(0)
    q(:,:,1)=qo_v;

    dq(:,:,1)=dqo_v;

    aq(:,:,1)=aqo;

%% set first u vector's point as u0.
u(:,:,1)=u0;


%% The iterative process
%Initializization
   dP_v=zeros(nn,1,length(t)); 
   dP_hat=zeros(nn,1,length(t));
   del_q=zeros(nn,1,length(t)); 
   del_dq=zeros(nn,1,length(t)); 
   del_aq=zeros(nn,1,length(t));  
%% Newmark's Method Iterative Process.
for i=1:length(t)-1

    dP_v(:,:,i)=P_v(:,:,i+1)-P_v(:,:,i); %Eqn 2.5.1.8
    dP_hat(:,:,i)=dP_v(:,:,i)+a*dq(:,:,i)+b*aq(:,:,i); %Eqn 2.5.1.9
    
    del_q(:,:,i)=Kh^(-1)*dP_hat(:,:,i); 
    %^Eqn 2.5.1.10
    del_dq(:,:,i)=(y/(beta*dt))*del_q(:,:,i)-(y/beta)*dq(:,:,i)+dt*(1-(y/(2*beta)))*aq(:,:,i);
    %^Eqn 2.5.1.11
    del_aq(:,:,i)=(1/(beta*dt^2))*del_q(:,:,i)-(1/(beta*dt))*dq(:,:,i)-(1/(2*beta))*aq(:,:,i);
    %^Eqn 2.5.1.12
    
    q(:,:,i+1)=q(:,:,i)+del_q(:,:,i); %Eqn 2.5.1.13
    dq(:,:,i+1)=dq(:,:,i)+del_dq(:,:,i); %Eqn 2.5.1.14
    aq(:,:,i+1)=aq(:,:,i)+del_aq(:,:,i); %Eqn 2.5.1.15
    
    u(:,:,i+1)=phi*q(:,:,i+1); %Eqn 2.5.1.16
end

%% Write Loop's results in easy-to-read table format. 

new_q=zeros(nn,length(t));
new_u=zeros(nn,length(t));
for i=1:length(t)
    for j=1:nn
    new_q(j,i)=q(j,1,i);
    new_u(j,i)=u(j,1,i);
    end
end

u=real(new_u');%sometimes it catches very small complex numbers.  
q=new_q';
% 
t=t';
Results_Table=[t,q,u]; %display t, qi(t),ui(t) column by column
%% Run FFT-Power Spectral Density 
if exist('FFT','var')
    figure(1)
    run('FFT_PSD') %do so on if loop
end

%% START Plotting everything. 

%% Plot #1 El_Centro Acceleration Time History
figure(2)
subplot(2,2,1) %% El_Centro Acceleration Time History northwest quadrant
plot(el_t,acc,'r'); %plot acceleration history  
hold on
gg=plot(xti,yi,'y.'); %plot the interpolated values of ground motion with yellow dots.
if Show_Time_Lapse<2
gh=plot(xti(1),yi(1),'bx','MarkerSize',10,'MarkerFaceColor','r','LineWidth',2); 
end
    %^^put a crosshair later to be animated at initial conditions
    grid
    xlabel('time (s)','FontSize',13)
    ylabel('Acceleration (g)','FontSize',13)
    title('El Centro Acceleration Time History','FontSize',16)
    xlim([t(1),t(end)])

%% Plot #2: Displacement by Newmark vs. Time
figure(2)
subplot(2,2,3) %% Displacement by Newmark vs. Time southwest quadrant
hold on

plot(t,u)
xlabel('t, [seconds]')
ylabel('u,[feet]')

%clever way of naming legends
for k=1:nn
    legendInfo{k}=['u' num2str(k)];
end
legend(legendInfo)
legend('Location','NorthEast')

%if user wants to see crosshair move along the graph...
if Show_Time_Lapse<2 
tq=plot(t(1),0,'bx','MarkerSize',15); %plots blue X in first subplot as time mechanism.
end

title('Displacement by Newmark-Avg. Acc','FontSize',15)
grid on
set(gcf,'Renderer','OpenGL');


% Plot #3 Chimney Plot
figure(2)
subplot(2,2,[2 4]) %each DOF's are plotted
hold on
index=0;
p=zeros(nn,1); 

%method of storing the dof values for them to be animated later.
for i=1:nn
    index=index+1;
p(index)=plot(u(1,i),h*i,'o','MarkerSize',15,'MarkerFaceColor','b'); %p are the blue dots.
%%^index a trick to hold onto plot values
end

set(gcf,'Renderer','OpenGL');

%window range is symmetric and fixed, so that the max displacmeent*110% is
%how wide the window will be on both ends. 

xmax=max(max(u));
xmin=min(min(u));
xwindow=max(xmax,abs(xmin));
xlim([-xwindow*1.1, xwindow*1.1]); %window range. 
ylim([0,nn*h*1.1]); %window range, start from base to top height*10% 

grid on
xlabel('u,displacement (ft)');
ylabel('Chimney height (ft)');
title('Chimney Displacement','FontSize',13);

%% Begin Animation
tic %measure how long it takes to run this loop
for i=1:length(t)
    for j=1:nn
set(p(j),'XData',u(i,j)); %set individual orb j to displace horizontally by u(i)
    end
    s=zeros(nn,1);
    s(1)=line([0,u(i,1)],[0,h],'Color','k','LineWidth',4);
    % ^draw chimney part from base to the 1st dof.

for d=2:nn
      s(d)=line([u(i,d-1),u(i,d)],[(d-1)*h,d*h],'Color','k','LineWidth',4);
      % ^draw chimney parts from dof d to dof d+1
end

   pause(dt) %animation runs at every dt frame. (makes animation run in real-time)
    if t(i)<t(end) %except for the very end of the loop
    delete(findobj(gca,'Color','black'));%find anything black ('chimney')
    %then delete it...way of animating chimney 
    end %ends the delete loop when program reaches its end. 
end
toc %measure how long it takes to run this loop.

