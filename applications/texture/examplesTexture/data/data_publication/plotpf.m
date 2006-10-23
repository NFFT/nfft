function plotpf(x, N1, r);

N2 = size(x, 1) / N1;
lines = floor(sqrt(N1));
cols = ceil(N1 / lines);

for i = 1:N1;
subplot(lines, cols, i);
rho = r(:, 1)';
theta = r(:, 2)';
c = x((1:N2) + (i-1)*N2)';

% Schmidt Projektion  - Flaechentreu
X = cos(rho) .* sqrt(2*(1-cos(theta)));
Y = sin(rho) .* sqrt(2*(1-cos(theta)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolierter plot:  -> Eingabe in Rechteckform
% X = [X;linspace(0, 1, N2)];
% Y = [Y;2*ones(1, N2)];
% Z = [c;ones(1, N2)];

% contour plot
% [CM,h] = contourf(X,Y,Z);
% set(h,'LineStyle','none')
% ColorMap(flipud(ColorMap('gray')));

% density plot
% pcolor(X,Y,Z);
% shading interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot as singular points 
diameter = 0.005;

hold on
cmin = min(reshape(c,[],1));
cmax = max(reshape(c,[],1));
c = 1+round((c-cmin) / (cmax-cmin) * 63);
cmap = colormap;

for i=1:numel(X)
	rectangle('Position',[X(i)-diameter/2,Y(i)-diameter/2,diameter,diameter],'Curvature',[1,1],...
		'FaceColor',cmap(c(i),:),'EdgeColor',cmap(c(i),:));
end

end;
