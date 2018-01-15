run('MDOF_ForcingFunction'); %user defined fun
%% MDOF Central Difference Method: Section 2.5.2 of Report
%Also mentioned in Chopra Table 15.2.1 pg 649
Tn=2*pi.*(1./(wn)); %natural angular frequency converted to natural period

w=nn; % #of modes to use. 
%% Stability Check
Tj=Tn(w); %set J'th period to the amount of modes user is using
if dt>0.1*Tj
    stability=0.1*Tj;%stability check 
    mystring='Soln unstable, lower time step dt below:';
    error(strcat(mystring,num2str(stability)));  
end
phi=phi(:,1:w);
K=K(1:w,1:w);
M=M(1:w,1:w);
C=C(1:w,1:w);


%% initial calculations
%1.1
    qo_v=zeros(1,length(nn)); 
    dqo_v=zeros(1,length(nn));
for i=1:w %for every column of phi vector; (for each mode shape)
    qo_v(i)=(phi(:,i)'*m*u0)/(phi(:,i)'*m*phi(:,i)); 
    %^ Eqn 2.5.2.1 get initial modal coordinate qo
    dqo_v(i)=(phi(:,i)'*m*v0)/(phi(:,i)'*m*phi(:,i)); 
    %^ Eqn 2.5.2.1 get initial modal "velocity" dqo
end
    qo_v=qo_v'; %transpose as column vector
    dqo_v=dqo_v'; %transpose as column vector
%1.2
Po=phi'*p_v(:,1,1); % Eqn 2.5.2.2 get initial P matrix normalized by phi.

%1.3: solve for double derivative of qo.
RHS=Po-C*dqo_v-K*qo_v; 
aqo=M^(-1)*RHS; % Eqn 2.5.2.3

%1.4 Select dt, and tr from User_Input_MDOF.

%1.5
q_1=qo_v-dt*dqo_v+(aqo*(dt)^2/2); % Eqn 2.5.2.4

%1.6 Get capital K hat
Kh=(1/dt^2)*M+(1/(2*dt))*C;%Eqn 2.5.2.5

%1.7 
a=(1/dt^2)*M(1:w,1:w)-(1/(2*dt))*C;  % Eqn 2.5.2.6
b=K-(2/(dt)^2)*M;  % Eqn 2.5.2.7

%%2.0 time step calculations
P_v=zeros(w,1,length(t));
for i=1:length(t)
    P_v(:,:,i)=phi'*p_v(:,:,i);  % Eqn 2.5.2.8
end
q=zeros(w,1,length(t));
q(:,:,1)=q_1;
dq(:,:,1)=dqo_v;
aq(:,:,1)=aqo;
u=zeros(nn,1,length(t));
u(:,:,1)=u0;
P_hat=zeros(w,1,length(t));

%% CDM iteration steps
for i=2:length(t)
    P_hat(:,:,i)=P_v(:,:,i)-a*q(:,:,i-1)-b*q(:,:,i);  % Eqn 2.5.2.9
    q(:,:,i+1)=Kh^(-1)*P_hat(:,:,i);  % Eqn 2.5.2.10
    u(:,:,i+1)=phi*q(:,:,i+1);  % Eqn 2.5.2.11
end

%% Write Loop's results in easy-to-read table format. 

new_q=zeros(w,length(t));
for i=1:length(t)
    for j=1:w
    new_q(j,i)=q(j,1,i);
    end
end
new_u=zeros(nn,length(t));
for i=1:length(t)
    for j=1:nn
    new_u(j,i)=u(j,1,i);
    end
end

u=real(new_u');%sometimes it catches very small complex numbers.  
q=new_q';
q=real(q);
% 
t=t';
Results_Table=[t,q,u]; %display t, qi(t),ui(t) column by column


%% Run FFT-Power Spectral Density if option chosen
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

title('Displacement by CDM','FontSize',15)
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

