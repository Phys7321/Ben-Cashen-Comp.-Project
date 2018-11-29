% Program to simulate protein hopping kinetics on DNA 
% GFP-LacI protein w/ DNA modeled as infinitely long, rigid cylinder w/
% constant radius

clear;

tic
delta = 0.267;      % Simulation step size (in angstroms)
tau = 4.46;         % Simulation step time (in picoseconds)
N = 1.1e+08;        % Number of steps per simulation (~(0.5 ms)/(4.46 ps))
R = 41.8;           % Capture radius (in angstroms)
N_s = 10000;        % Number of simulations
N_accept = 0;
N_reject = 0;

% Allocating space for measurement arrays
z = zeros(N+1,1);
r = zeros(2,N+1);
r_mag = zeros(N+1,1);
dist = zeros(N_s,1);
n = zeros(N_s,1);
time = zeros(N_s,1);
h = zeros(N_s,1);
hop_num = zeros(N_s,1);

rng('shuffle')

for j = 1:N_s
    
% Initial position (in angstroms)
x = 42.1;
y = 0;
z(1) = 0;
r(:,1) = [42.1; 0];
r_mag(1) = 42.1;

% Drawing step values from Gaussian distribution w/ stdev of delta and mean
% of zero (column vectors of size N)
dx = normrnd(0,delta,[N,1]);   
dy = normrnd(0,delta,[N,1]);
dz = normrnd(0,delta,[N,1]);

% Letting the protein experience that sweet feeling of unbiased diffusion
for i = 1:N
    x = x + dx(i);
    y = y + dy(i);
    z(i+1) = z(i) + dz(i);
    r(:,i+1) = [x; y];
    r_mag(i+1) = norm(r(:,i+1));
    d = r(:,i+1) - r(:,i);                          % Displacement vector between consecutive iterations
    dot1 = d(1)*r(1,i) + d(2)*r(2,i);               % dot(d,r(:,i))
    dot2 = d(1)*r(1,i+1) + d(2)*r(2,i+1);           % dot(d,r(:,i+1))
    cross = d(1)*r(2,i) - d(2)*r(1,i);              % cross(d,r(:,i))
    
    if (dot1/dot2 <= 0 && norm(cross)/norm(d) <= R) || (r_mag(i+1) && r_mag(i) <= R)
        break
    end
end

if i < N                                            % Only accept if association occurs within set time frame
    N_accept = N_accept + 1;                        % Number of hops 
    h(N_accept) = max(r_mag(1:i+1)) - 41.8;         % Max. hopping height
    dist(N_accept) = abs((z(i)+z(i+1))/2);          % Hopping distance along DNA
    n(N_accept) = i;                                % Number of steps per hop
    time(N_accept) = n(N_accept)*tau*10^(-6);       % Hopping time (microseconds)
else
    N_reject = N_reject + 1;                        % Dissociation events
    hop_num(N_reject) = N_accept;                   % Number of hops per diffusion trajectory
end
end

hop_num(N_reject+1) = N_accept - hop_num(N_reject);
hop_num(2:N_reject) = diff(hop_num(1:N_reject));

% Creating bin edge vector for histogram
bin1 = logspace(-3,3,30);
bin2 = logspace(-1,3,30);
bin3 = logspace(-6,3,45);
bin4 = linspace(0,900,30);

% Average values after N_s hopping simulations
avg_h = sum(h)/N_accept;                            % Avg. max hopping height
avg_dist = sum(dist)/N_accept;                      % Avg. hopping distance
avg_n = sum(n)/N_accept;                            % Avg. steps per hop
avg_time = sum(time)/N_accept;                      % Avg. time per hop
avg_hop_num = sum(hop_num)/(N_reject);              % Avg. hops per trajectory

toc

% Distribution of hopping distance
figure(1)
histogram(dist(1:N_accept),bin1)
set(gca,'XScale','log')
set(gca,'YScale','log')
title('Distribution of hopping distance (GFP)')
xlabel('Hopping distance (\AA)')
ylabel('Count')

% Distribution of maximum hopping height
figure(2)
histogram(h(1:N_accept),bin2)
set(gca,'XScale','log')
set(gca,'YScale','log')
title('Distribution of max. hopping heights (GFP)')
xlabel('Hopping height (\AA)')
ylabel('Count')

% Distribution of hopping time
figure(3)
histogram(time(1:N_accept),bin3)
set(gca,'XScale','log')
set(gca,'YScale','log')
title('Distribution of hopping time (GFP)')
xlabel('Hopping time ($\mu$s)')
ylabel('Count')

% Distribution of hops per diffusion trajectory
figure(4)
histogram(hop_num(1:N_reject+1),bin4);
title('Distribution of hops per diffusion trajectory (GFP)')
xlabel('Hops per trajectory')
ylabel('Count')