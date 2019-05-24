function [ output ] = fun_clus_noRW( kact, bzf)

%% Variables:
% kact = basal methylation rate
% bzf = ZF specificity multiplier

% % Test Variables
% kact = 5*10^-4;
% bzf = 10;

%% Rate constants (p = parameter structure array)

p.kdecay = 0.05;                            % Mark Turnover Rate

% Initiator properties
p.Init_kmark = bzf*kact;                    % ZF-specific Initiation
p.Init_basal = kact;                        % Basal methylation by Initiator


%% Initial state
numSites = 63;                              % number of GATC sites
tspan = [0, 5000];                          % length of simulation (minutes)
x0 = zeros(1,numSites);                     % initial state of lattice


%% Specify reaction network
pfun = @lattice_clus_noRW;

stoich_matrix = [ 1  zeros(1,numSites-1);   % nucleation:    U --  kmark --> M
                  eye(numSites)*-1          % turnover:      M -- kdecay --> U
                  eye(numSites); ];         % basal marking: U --  kact --> M

              
%% Run Simulation
numSim = 1000;                              % Number of simulations

% Initialize Storage vectors
stor_end = zeros(numSim,numSites);          % store steady state value

% Run simulations
parfor i=1:numSim
    [t,x] = directMethod(stoich_matrix, pfun, tspan, x0, p);
    stor_end(i,:) = x(end,:);               % Store end of simulation
end


%% Heatmap plot of modification profile

L = numSim;   % Number of simulations to plot

figure(1)
imagesc(stor_end(1:L,1:63))
map = [ 256 256 256;            % Colormap (white -> black)
        0   0   0]/256;
colormap(map)
xlim([1 numSites])
set(gca,'FontSize',18)
xlabel('GATC position');
ylabel('steady state profile over 1000 sims'); 


%% Mean distribution of modification along array
x = 1:numSites;

figure(2)
boundedline(x, mean(stor_end), std(stor_end).^2,'k-'); hold on
xlim([1 numSites])
ylim([0 1.1])
set(gca,'FontSize',22)
xlabel('GATC position')
ylabel('mean mark density')
pbaspect([1.2 1 1])


%% Output Parameters
output.stor_end = stor_end;     % Final methylation profiles at end of each simulation
output.params = p;              % Store Parameters

end

