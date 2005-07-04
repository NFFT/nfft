it_max = 10;
d_phi = 720
d_theta = 360
for k=0:it_max
  filename = sprintf('it%02d.dat',k);
  f = fopen(filename,'r');
  topo = fscanf(f,'%E\n',[d_theta,d_phi]);
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
  filename = sprintf('it%02d.jpg',k);
  imwrite(toporgb,filename,'jpg','Quality',100);  
  %F(k+1) = getframe(gcf)
  %imwrite(A,map,filename,'tga','XResolution',720,'YResolution',360)  
end
%movie(F)
  filename = sprintf('it_orig.dat',k);
  f = fopen(filename,'r');
  topo = fscanf(f,'%E\n',[d_theta,d_phi]);
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
  filename = sprintf('it_orig.jpg',k);
  imwrite(toporgb,filename,'jpg','Quality',100);  
  %F(k+1) = getframe(gcf)
  %imwrite(A,map,filename,'tga','XResolution',720,'YResolution',360)  
