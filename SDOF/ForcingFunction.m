%% SAMPLE FORCING FUNCTION INPUT #1: Subject user-created structure to a
%hypothetical sinusoidal forcing function

%% User Input
ptr=0.6; %forcing function time range [s]
po=10; %forcing function amplitude [kips]
wf=pi/0.6; %forcing function frequency [rad/s]
%% END User Input

p=zeros(1,length(t)); %initialize p vector,
for i=1:length(t) 
    if t(i)<=ptr %during the interval the forcing function is applied...
        p(i)=po*sin(wf*t(i)); % the forcing function is this...
        % ^^same forcing function from Chopra example 5.1 (206)
    else
        p(i)=0; %outside of the applied force range, the force=zero.
    end
end


%% SAMPLE FORCING FUNCTION INPUT #2: Subject user-created structure to the 
% %El Centro Earthquake

%%BEGIN User Input
% run('El_Centro_Ground_Motion') %feed Data downloaded from NCEES website.
% %%END User Input
% 
% b=a'; %transpose the ground history data...  
% zz=b(:); %so that the data turns into neat column vector;
% zz(length(zz))=[]; %omit last value in data;0 added to square the matrix
% acc=[0;zz]; %@ t=0, acceleration is 0, motion data starts at t=0.02; add 0
% el_t=0:0.02:31.18; %time vector goes from 0 to 31.18s with time step .02s;
% x=el_t'; %El Centro Time vector converted to column vector
% y=acc; %El Centro acc vector is already a column vector.
% xi=t'; %convert user input time vector into column vector
% yi=interp1q(x,y,xi); %linearly interpol. to attain acc pts spaced dt apart.  
% yi(isnan(yi))=0; %set values after final data in ground motion acc=0.  
% g=32.2*12; %gravity in inches/s^2
% ag=yi*g; %ground motion acc. a(t);
% p=-m*ag; %external force due to ground motion in SDOF system  
% p=p'; %revert p vector into row vector to match the code consistency. 



    

