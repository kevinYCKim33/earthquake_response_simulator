% FFT=1; %Simple way to run up to FFT only.  
% run('MDOF_Run');

fs=1/dt; %sample frequency;
%t=already have time vector--a column vector at this point.

catch_u=u(:,1); %feed dof1's u into the fft function.
X=fft(catch_u);
%a complex number at this point
X_mag=abs(X); %gets the norm of real and i function to get pure pos #'s

figure(1)
freq=linspace(0,fs,length(X_mag)); %set x-scale normalized to frequency
plot(freq,X_mag/(.5*length(X_mag))); %normalize where peaks are displayed in magnitude.
title('Power Spectral Density','FontSize',15)
xlabel('Frequency Hz','FontSize',13)
ylabel('ft','FontSize',13)
grid on
fNy=.5*fs;
fn=wn/(2*pi);
xlim([0,abs(fn(nn))*1.1]); %set the x-axis window as being slightly
%longer than the last frequency value
hold on
for i=1:w %all the modal frequencies...
   ef_n=plot(fn(i),0,'rx','MarkerSize',15); %are plotted with an 'x' on FFT window
end
legend(ef_n,'theoretic mode shape frequencies','FontSize',13)
hold off
% dash_line=line([t(1),t(1)],get(gca,'YLim'),'Color','black','LineStyle','--'); 
% put a dashed line to be moving across as time goes 
