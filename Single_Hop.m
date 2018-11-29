% Program to simulate and animate protein diffusion around DNA 
% Capture criteria removed 

clear;

delta = 0.5;                % Simulation step size made artificially high for visualization (in angstroms)
tau = 4.46;                 % Simulation step time (in picoseconds)
N = 300;                    % Number of steps per simulation
R = 41.8;                   % Capture radius (in angstroms)

% Allocating space for measurement arrays
x = zeros(N+1,1);
y = zeros(N+1,1);
z = zeros(N+1,1);
F(N+1) = struct('cdata',[],'colormap',[]);

% Defining periodic radius of cylinder (Z-DNA)
t = linspace(0, 2*pi);
cyl = [R*cos(t); R*sin(t)];

rng('shuffle')

% Initial position (in angstroms)
x(1) = 42.5; 
y(1) = 0;
z(1) = 0;

% Drawing step values from Gaussian distribution w/ stdev of delta and mean
% of zero (column vectors of size N)
dx = normrnd(0,delta,[N,1]);   
dy = normrnd(0,delta,[N,1]);
dz = normrnd(0,delta,[N,1]);

% Letting the protein experience that sweet feeling of unbiased diffusion
vidfile = VideoWriter('testmovie.mp4','MPEG-4');
open(vidfile);
for i = 1:N
    x(i+1) = x(i) + dx(i);
    y(i+1) = y(i) + dy(i);
    z(i+1) = z(i) + dz(i);
    surf([cyl(1,:); cyl(1,:)], [cyl(2,:); cyl(2,:)], [-100*ones(1,size(cyl,2)); 100*ones(1,size(cyl,2))],'FaceColor','w')
    hold on
    c = z;
    scatter3(x,y,z,3,c)
    view(20,20)
    %axis equal
    drawnow limitrate
    F(i) = getframe(gcf);
    writeVideo(vidfile,F(i));
    hold off
    rotate3d on
end
close(vidfile)
