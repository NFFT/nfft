N=256;
border_eps=0.2;

M=8345;

load vol87.dat -ascii
input_data=vol87(2:end,:);

x_range=max(input_data(:,1))-min(input_data(:,1));
input_data(:,1)=(input_data(:,1)-min(input_data(:,1))) /x_range*(1-2*border_eps)-0.5+border_eps;

y_range=max(input_data(:,2))-min(input_data(:,2));
input_data(:,2)=(input_data(:,2)-min(input_data(:,2))) /y_range*(1-2*border_eps)-0.5+border_eps;

save input_data.dat -ascii -double -tabs input_data

!make

!glacier 256 8345 > output_data.dat

load output_data.dat
f_hat=output_data(:,1)+i*output_data(:,2);

figure(1);
C=ifftshift(real(fft2(fftshift(reshape(f_hat,N,N)))));
x=-0.5:1/N:0.5-1/N;
[X,Y]=meshgrid(x,x);
surfl(X,Y,C);
axis([-0.5+border_eps,0.5-border_eps,-0.5+border_eps,0.5-border_eps,min(input_data(:,3)),max(C(:))]);
colormap(bone);
shading interp
axis off
view([30,30])
print glacier1.eps -deps

figure(2);
contour(X,Y,C,sort(input_data(:,3)));                                                    
axis([-0.5+border_eps,0.5-border_eps,-0.5+border_eps,0.5-border_eps]);
axis off
hold on
plot(input_data(:,1),input_data(:,2),'k.');
colormap(gray(256));
hold off
print glacier2.eps -deps
