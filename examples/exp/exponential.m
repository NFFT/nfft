k=2.5;
x=0.5;

system(['./exp  '  num2str(k) ' ' num2str(x) ' > out.dat']);

load out.dat

hold on;
plot(out(:,1),out(:,2),'r-');
plot(out(:,1),out(:,3),'b-');
plot(out(:,1),out(:,4),'g-');
