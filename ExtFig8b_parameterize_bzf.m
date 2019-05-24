clc, clear all, close all

%% ANALYZE MEASURED DATA

% Load data
load Dam_properties.mat

% Get mean values
spec_mean = mean(spec')';       % Mean ZF-dam/mch-dam frac methylated
mch_dam = mean(act')';          % Mean mch-dam frac methylated
ZF_dam = mch_dam.*spec_mean;    % Mean ZF-dam frac methylated


%% FUNCTIONS TO PREDICT bzf (ZF-specificity multiplier)

% Parameters
k_dec = 0.05;                           % Mark Turnover rate (scales with doubling time)
k_act_range = logspace(-6,1,100);       % Range of k_act rates

% Fraction Methylated (probability bound) - equations in supp text
pbound = @(k_act,k_dec) k_act./(k_dec+k_act);                                     % frac meth for mch-Dam
pbound_withZF = @(bzf,k_act,k_dec) (bzf.*k_act+k_act)./(k_dec+bzf.*k_act+k_act);  % frac meth for ZF-Dam


%% PLOT

figure
    loglog(pbound(k_act_range,k_dec),pbound_withZF(1,k_act_range,k_dec),'r--'); hold on   
    loglog(pbound(k_act_range,k_dec),pbound_withZF(10,k_act_range,k_dec),'k--'); hold on   
    loglog(pbound(k_act_range,k_dec),pbound_withZF(100,k_act_range,k_dec),'b--'); hold on  
    
    loglog(mch_dam,ZF_dam,'.','MarkerSize',30); hold on
    
    loglog([10^-4 1.1],[10^-4 1.1],'k--')
    set(gca,'FontSize',18)
    xlabel('mch-Dam')
    ylabel('ZF-Dam')
    legend('bzf=1','bzf=10','bzf=100','Location','SouthEast')
    ylim([10^-4 1.1])
    xlim([10^-4 1.1])
    pbaspect([1.1 1 1])
    