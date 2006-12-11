N=256;
border_eps=0.1;

M=8345;

load vol87.dat -ascii
input_data=vol87(2:end,:);

x_range=max(input_data(:,1))-min(input_data(:,1));
input_data(:,1)=(input_data(:,1)-min(input_data(:,1))) /x_range*(1-2*border_eps)-0.5+border_eps;

y_range=max(input_data(:,2))-min(input_data(:,2));
input_data(:,2)=(input_data(:,2)-min(input_data(:,2))) /y_range*(1-2*border_eps)-0.5+border_eps;

% Resort samples randomly
P=randperm(M);
input_data=input_data(P,:);

M_cv_start=200;
M_cv_step=200;
M_cv_end=1000;

save input_data.dat -ascii -double -tabs input_data

system(sprintf('./glacier %d %d %d %d %d> output_data_cv.tex',N,M,M_cv_start,M_cv_step,M_cv_end));
