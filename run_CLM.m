% Simulates Charge Lattice Model and Associated Lattice Distortions
% 2D charge lattice is thermally disordered using Monte Carlos method
% using Metropolis-Hastings alogrithm with Lennard-Jones interactions.
% Associated Lattice Distortions is calculated by displacing crystal 
% lattice points along the nearest charge lattice sites.

% This script produces two end result variables
% 1. pos_Chas: numSiteCha x 2 x numT
%    Charge lattice point coordinates at each temperature
% 2. pos_Cry_PLD: numSiteCry x 2 x numT
%    Crystal lattice coorindates with associated lattice distortion 

% Simulated results are saved in results/seed_number.mat


%% Simulation Setup.
boolDraw = true; % Set to true to plot simulation progress. Can be very slow.
mul = 3; % Number of x,y tiling
numIter = 10^3; % Number of Monte Carlo Iterations for Temperature

Tmax = 1000; % Max temperature in unit of V_min. V_min is minimum energy of LJ potential
numT = 15;  % Number of Temperature Steps to simulate between (0, Tmax);

%% Avoid modifying anything below here

%% Shuffle random seed
rng('shuffle')
sd = rng; 

%% Tools
getR = @(x) sqrt( x(:,1).^2 + x(:,2).^2 );

%% Lattice Grid Setup
% Cry stands for crystal lattice
% Cha stands for charge superlattice

% Lattice Vectors
a_Cry = 3; th_Cry = 0;
a1_Cry = a_Cry*[cosd(th_Cry), sind(th_Cry)]; 
a2_Cry = a_Cry*[cosd(th_Cry+120),sind(th_Cry+120)];

a_Cha = sqrt(13)*a_Cry; th_Cha = th_Cry;
a1_Cha = a_Cha*[cosd(th_Cha), sind(th_Cha)]; 
a2_Cha = a_Cha*[cosd(th_Cha+120), sind(th_Cha+120)];

% Populate Big lattice and cut off later
[m1,m2] = meshgrid( -800:800 );
[n1,n2] = meshgrid( -400:400 );

m1 = m1(:); m2 = m2(:);
n1 = n1(:); n2 = n2(:);

pos_Cry = m1*a1_Cry + m2*a2_Cry;
pos_Cha = n1*a1_Cha + n2*a2_Cha;

% Periodic Boundary size, Roughly square and commensurate with lattice
xCutoffCha = a_Cha*21*mul; yCutoffCha = xCutoffCha/7*4*sqrt(3);

% Round off values to avoid precision error
rd = 4;
a_Cry = round(a_Cry,rd); a1_Cry = round(a1_Cry, rd); a2_Cry = round(a2_Cry,rd);
a_Cha = round(a_Cha,rd); a1_Cha = round(a1_Cha, rd); a2_Cha = round(a2_Cha, rd);
pos_Cry = round(pos_Cry,rd); pos_Cha = round( pos_Cha, rd);
xCutoffCha = round( xCutoffCha, rd); yCutoffCha = round( yCutoffCha, rd);

% Remove Charge Lattice Points outside of periodic boundary
pos_Cha( pos_Cha(:,1) >=  xCutoffCha,: ) = [];
pos_Cha( pos_Cha(:,1) <   0         ,: ) = [];
pos_Cha( pos_Cha(:,2) >=  yCutoffCha,: ) = [];
pos_Cha( pos_Cha(:,2) <   0         ,: ) = [];

% Remove Crystal Lattice Points outside of cutoff radius from center of
% periodic boundary
rpos_Cry = sqrt( (pos_Cry(:,1)-xCutoffCha/2).^2 + (pos_Cry(:,2)-yCutoffCha/2).^2);
pos_Cry( rpos_Cry > xCutoffCha/2*0.7, : ) = [];

numSiteCry = length(pos_Cry);
numSiteCha = length(pos_Cha);

%% Monte Carlo Setup

% Truncated LJ Potential
sgLJ = 1.5*a_Cha*(2^(-1/6)); ep = 5;
A = 4*ep*sgLJ^12;
B = 4*ep*sgLJ^6;
r_c = sqrt(3)*a_Cha;

% Iteration Setup
kBTs = linspace(0, Tmax, numT);

iter = numIter*ones([numT,1]);
iter(1) = 0;

dr_max = a_Cha*0.1; % Maximum Displacement allowed per iteration

pos_Cry_PLD = zeros([size(pos_Cry),numT]);
pos_Chas    = zeros([size(pos_Cha),numT]);
pldScal = 0.1;
%% Iterate

if boolDraw
    figure
    hold on
    plot([0, xCutoffCha, xCutoffCha,0,0],[0,0,yCutoffCha,yCutoffCha,0],'k-')
    sh = scatter(pos_Cha(:,1), pos_Cha(:,2),'k.');
    axis equal
    xlim(xCutoffCha*[-0.05, 1.05]);
    ylim(yCutoffCha*[-0.05, 1.05]);
    drawnow
end

for indT = 1:numT
    kBT = kBTs(indT);
    tic
    fprintf('%d/%d\n',indT,numT)
    
    for ct = 1:iter(indT)
        
        for indCha = 1:numSiteCha
            pos_cur = pos_Cha(indCha,:);
            
            dPosR = dr_max*(2*rand-1);
            dPosT = 2*pi*rand;
            dPos = round(dPosR*[cos(dPosT),sin(dPosT)],rd);
            %dPos = round(dr_max*(2*rand(1,2)-1),rd);
            
            dR_nns = pos_cur - pos_Cha;
            dR_nns(:,1) = dR_nns(:,1) - round(dR_nns(:,1) / xCutoffCha) * xCutoffCha;
            dR_nns(:,2) = dR_nns(:,2) - round(dR_nns(:,2) / yCutoffCha) * yCutoffCha;
            
            r_nns = getR(dR_nns);
            dR_nns( r_nns > r_c | r_nns <= 0, :) = [];
            
            
            dR_nns_new = dR_nns + dPos;
            r_nns_cur = getR(dR_nns);
            r_nns_new = getR(dR_nns_new);

            dE = sum(V_LJ(r_nns_new, A, B,r_c)) - sum(V_LJ(r_nns_cur, A, B,r_c));

            if rand() < exp(-dE/(kBT))
                pos_new = pos_cur + dPos;
                pos_Cha(indCha,1) = pos_new(1) - floor(pos_new(1)/xCutoffCha)* xCutoffCha;
                pos_Cha(indCha,2) = pos_new(2) - floor(pos_new(2)/yCutoffCha)* yCutoffCha;
            end
        end
        if boolDraw
            sh.XData = pos_Cha(:,1); sh.YData = pos_Cha(:,2);
            drawnow
        end
    end
    pos_Chas(:,:,indT)    = pos_Cha;
    pos_Cry_PLD(:,:,indT) = calcPLD_gaussian(pos_Cry, pos_Cha, a_Cha/2, a_Cha/2, pldScal);
    
    toc
end

%% Save Data
if ~isfolder('results')
    mkdir('results')
end
save( sprintf('results/%d.mat',(sd.Seed)) )