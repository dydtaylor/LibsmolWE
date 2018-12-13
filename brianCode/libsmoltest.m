



tic
load trajectories.mat
for i = 1:10
tr1 = trajectories(1,:);
%tr = reshape(trajectories,[size(trajectories,2)*size(trajectories,1) 1]);
nc = size(tr1,2);
a=evacuation(tr1,nc)
toc

end
