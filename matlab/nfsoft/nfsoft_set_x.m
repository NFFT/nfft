%Set rotations in Euler angles for nfsoft
% x(1,:) = alpha
% x(2,:) = beta
% x(3,:) = gamma
function nfsoft_set_x(plan,x)
nfsoftmex('set_x',plan,x)