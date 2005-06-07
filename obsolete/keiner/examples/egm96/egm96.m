k_min = 0
k_max = 7
d_phi = 720
d_theta = 360
for k=k_min:k_max
  filename = sprintf('egm96%05d.dat',k);
  f = fopen(filename,'r');
  topo = fscanf(f,'%f\n',[d_theta,d_phi]);
  imagesc(topo)
  fclose(f);
  F(k+1) = getframe(gcf)
end
movie(F)
