function [ output ] = fun_clus_withRW( kact_synI, bzf, kact_RW, bdpn1)

%% Variables:
% kact_SynI = basal methylation rate of SynI
% bzf = ZF specificity multiplier
% kact_RW = basal methylation rate of SynRW
% bzf = DpnI specificity multiplier

% % Test Variables
% kact_synI = 5*10^-4;
% bzf = 10;
% kact_RW = 5*10^-4.2;
% bdpn1 = 100;

%% Rate constants (p = parameter structure array)

p.kdecay = 0.05;                            % Mark Turnover rate

% Initiator properties
p.Init_kmark = bzf*kact_synI;               % ZF-specific Initiation by SynI
p.Init_basal = kact_synI;                   % Basal methylation by SynI
p.Init_prop  = 1*10^-3;                     % Propogation by SynI

% Read-Writer properties
p.RW_basal = kact_RW;                       % Basal methylation by SynRW
p.RW_prop  = bdpn1*kact_RW;                 % Propogation rate by SynRW

%% Initial state
numSites = 63;                              % number of GATC sites
tspan = [0, 5000];                          % length of simulation (minutes)
x0 = zeros(1,numSites);                     % initial state of lattice

%% Specify reaction network
pfun = @lattice_clus_withRW;

stoich_matrix = [ 1  zeros(1,numSites-1);   % nucleation:  U --  Init_kmark --> M
                  eye(numSites);            % propogation: U -- RW_prop --> M
                  eye(numSites)*-1          % turnover:    M -- kdecay --> U
                  eye(numSites); ];         % basal:       U --  kact_synI + kact_RW --> M

%% Run Simulation
numSim = 1000;                              % Number of simulations

% Initialize Storage vectors
stor_end = zeros(numSim,numSites);          % store steady state value
stor_n15 = zeros(numSim,1);                 % store mean number of marks

% Run simulations
parfor i=1:numSim
    [t,x] = directMethod(stoich_matrix, pfun, tspan, x0, p);
    stor_end(i,:) = x(end,:);                       % Store end of simulation
    stor_n15(i) = mean(sum(x(:,end-14:end)'));      % Store 15 marks farthest from nucleation site
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
output.stor_end = stor_end;         % Final methylation profiles at end of each simulation
output.params = p;                  % Store Parameters
output.md15 = mean(stor_n15);       % Store mean mark density of 15 distant GATC sites
output.std15 = std(stor_n15);       % Store std of mark density of 15 distant GATC sites

end

