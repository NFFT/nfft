N=20
rec=16;
M=10

!make

h=1.2/N;
X=my_rand(M,h);
f=[0*X(X<=0);0.5^-0.5*X(X>0).^0.5]; 	%rand(M,1);
input_data=[X,f];
save input_data.dat -ascii -double -tabs input_data

!interpolation_1d 20 10 0 > output_data.dat
load output_data.dat
f0_hat=output_data(:,1)+i*output_data(:,2);

!interpolation_1d 20 10 1 > output_data.dat
load output_data.dat
f1_hat=output_data(:,1)+i*output_data(:,2);

!interpolation_1d 20 10 2 > output_data.dat
load output_data.dat
f2_hat=output_data(:,1)+i*output_data(:,2);

!interpolation_1d 20 10 3 > output_data.dat
load output_data.dat
f3_hat=output_data(:,1)+i*output_data(:,2);

!interpolation_1d 20 10 4 > output_data.dat
load output_data.dat
f4_hat=output_data(:,1)+i*output_data(:,2);

!interpolation_1d 20 10 5 > output_data.dat
load output_data.dat
f5_hat=output_data(:,1)+i*output_data(:,2);


N_rec=rec*N;
X_rec=((-0.5):1/N_rec:(0.5-1/N_rec))';
f_known=[0*X_rec(X_rec<=0);0.5^-0.5*X_rec(X_rec>0).^0.5];
f0_rec=ifftshift(real(fft(fftshift([zeros((N_rec-N)/2,1);f0_hat;zeros((N_rec-N)/2,1)]))));
f1_rec=ifftshift(real(fft(fftshift([zeros((N_rec-N)/2,1);f1_hat;zeros((N_rec-N)/2,1)]))));
f2_rec=ifftshift(real(fft(fftshift([zeros((N_rec-N)/2,1);f2_hat;zeros((N_rec-N)/2,1)]))));
f3_rec=ifftshift(real(fft(fftshift([zeros((N_rec-N)/2,1);f3_hat;zeros((N_rec-N)/2,1)]))));
f4_rec=ifftshift(real(fft(fftshift([zeros((N_rec-N)/2,1);f4_hat;zeros((N_rec-N)/2,1)]))));
f5_rec=ifftshift(real(fft(fftshift([zeros((N_rec-N)/2,1);f5_hat;zeros((N_rec-N)/2,1)]))));

figure(1); h=plot(X,f,'ko',X_rec,f0_rec,'k',X_rec,f_known,'k--'); axis([-0.5,0.5,-0.1,1.1]);
set(h,'LineWidth',2.0); set(h,'Markersize',10); set(gca,'FontSize',22); print interpolation_1d_0.eps -deps
figure(2); h=plot(X,f,'ko',X_rec,f1_rec,'k',X_rec,f_known,'k--'); axis([-0.5,0.5,-0.1,1.1]);
set(h,'LineWidth',2.0); set(h,'Markersize',10); set(gca,'FontSize',22); print interpolation_1d_1.eps -deps
figure(3); h=plot(X,f,'ko',X_rec,f2_rec,'k',X_rec,f_known,'k--'); axis([-0.5,0.5,-0.1,1.1]);
set(h,'LineWidth',2.0); set(h,'Markersize',10); set(gca,'FontSize',22); print interpolation_1d_2.eps -deps
figure(4); h=plot(X,f,'ko',X_rec,f3_rec,'k',X_rec,f_known,'k--'); axis([-0.5,0.5,-0.1,1.1]);
set(h,'LineWidth',2.0); set(h,'Markersize',10); set(gca,'FontSize',22); print interpolation_1d_3.eps -deps
figure(5); h=plot(X,f,'ko',X_rec,f4_rec,'k',X_rec,f_known,'k--'); axis([-0.5,0.5,-0.1,1.1]);
set(h,'LineWidth',2.0); set(h,'Markersize',10); set(gca,'FontSize',22); print interpolation_1d_4.eps -deps
figure(6); h=plot(X,f,'ko',X_rec,f5_rec,'k',X_rec,f_known,'k--'); axis([-0.5,0.5,-0.1,1.1]);
set(h,'LineWidth',2.0); set(h,'Markersize',10); set(gca,'FontSize',22); print interpolation_1d_5.eps -deps

