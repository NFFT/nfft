k_min = 0
k_max = 357
d_phi = 720
d_theta = 360
for k=k_min:k_max
  k
  filename = sprintf('egm96%05d.dat',k);
  f = fopen(filename,'r');
  topo = fscanf(f,'%f\n',[d_theta,d_phi]);
  fclose(f);
  ma = max(max(topo));
  mi = min(min(topo));
  topo = uint8((64/ma)*(topo-mi));
  %ma = max(max(topo))
  %mi = min(min(topo))
  map = colormap;
  %whos topo
  %whos map
  toporgb = ind2rgb(topo,map);
  filename = sprintf('egm96%05d.jpg',k);
  imwrite(toporgb,filename,'jpg','Quality',100);  
  %F(k+1) = getframe(gcf)
  %imwrite(A,map,filename,'tga','XResolution',720,'YResolution',360)  
end
%movie(F)
