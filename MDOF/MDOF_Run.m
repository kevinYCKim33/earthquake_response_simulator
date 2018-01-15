
clc
clear all;
close all;
run('User_Input_MDOF')
nn=length(m); %number of nodes
t=(0:dt:tr); %time vector [s]
%% create mass matrix
m_temp=zeros(nn,nn);
for i=1:nn
    m_temp(i,i)=m(i); %put m(i) along trace of square matrix  
end
m=m_temp;  %m matrix is now a diagonal matrix.

run('MDOF_ForcingFunction')
syms wsq %w^2, the effective eignvalue of modal frequencies.

%% create statically condensed k matrix

%% ktt
ktt=zeros(nn,nn);
a=24*EI/(h^3);
b=-12*EI/(h^3);
c=12*EI/(h^3);
for i=1:nn-1
    ktt(i,i)=a;
end
ktt(nn,nn)=c;
for i=2:nn
    ktt(i,i-1)=b;
end

for i=1:nn-1
    ktt(i,i+1)=b;
end

%% kot
kot=zeros(nn,nn);
d=6*EI/(h^2);
kot(nn,nn)=d;
for i=2:nn
    kot(i,i-1)=-d;
end

for i=1:nn-1
    kot(i,i+1)=d;
end

%% kto
kto=zeros(nn,nn);

kto(nn,nn)=d;
for i=2:nn
    kto(i,i-1)=d;
end

for i=1:nn-1
    kto(i,i+1)=-d;
end

%% koo
% %bottom-right quadrant
koo=zeros(nn,nn);

e=8*EI/h;
f=2*EI/h;
g=4*EI/h;

for i=1:nn-1
    koo(i,i)=e;
end
koo(nn,nn)=g;
for i=2:nn
    koo(i,i-1)=f;
end

for i=1:nn-1
    koo(i,i+1)=f;
end

ktt_hat=ktt-kot.'*koo^(-1)*kot; %statically condensed ktt^
k=ktt_hat; %for convenience, just call ktt^ k from now on. 


%% Get Modal Freqencies by Eigenvalue process
A=k-wsq*m; 

detA=det(A);

eigs=solve(detA==0,wsq); %solve for eigenvalues--won't work on R2010a
eigs=double(eigs); %get it into numbers and decimals
wn=sqrt(eigs); %take square-root to get modal frequencies
wn=sort(wn,'ascend'); %arrange modal frequencies from smallest to largest

%% Find Mode Shapes
phi=zeros(nn,nn);

for kk=1:nn 
    %
    B_temp=k-wn(kk)^2*m; 
    C_temp=B_temp(2:nn,2:nn); 
    D_temp=-B_temp(2:nn,1);
    f_temp=C_temp\D_temp;
    phi_temp=[1;f_temp]; %set the first entry as 1;
    for ii=1:nn
        phi(ii,kk)=phi_temp(ii);%store phi as a 2D matrix.
    end
end

%% Damping Matrix by Rayleigh Damping (Eqn. 11.4.10 from Chopra p493)
rao=(xi*2*wn(1)*wn(2))/(wn(1)+wn(2)); %Rayleigh coefficient a0.
ra1=(xi*2)/(wn(1)+wn(2)); %Rayleigh coefficient a1;
lil_c=rao*m+ra1*k; %damping matrix (11.4.7)  


%% uncouple m,k,c matrices by taking the triple product.
M=phi.'*m*phi;
K=phi.'*k*phi;
C=phi.'*lil_c*phi;
% P=phi.'*p_v(:,1,2); %  NOT CODE AT THIS POINT!


%% MDOF Analysis Type chosen from here.  
switch lower(Analysis_Type) %turns the user-inputs into lowercase letters

case 'cdm' %for Central Difference Method Case...
    if Ground_Motion>1 %if user is using smooth forcing function
        run('MDOF_Central_Difference_Method') %run this program
    else
        run('MDOF_Central_Difference_Method_Ground_Motion') %run this for 
        %running El_Centro Ground Motion
    end
case 'newmark' %for Newmark's Method Case
    if Ground_Motion>1%if user is using smooth forcing function
        run('MDOF_Newmark') %% run this program (pg.650; Table 15.2.2)
    else
        run('MDOF_Newmark_Ground_Motion')%run this for 
        %running El_Centro Ground Motion
    end
   
end



