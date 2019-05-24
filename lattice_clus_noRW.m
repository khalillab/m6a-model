function a = lattice_clus_noRW(x, p)
% Return reaction propensities given current state x
% x = current state of lattice
% p = structure array with rate parameters

%% nucleation
start_site = x(1);
if start_site == 0
    nuc_bool = 1;
else
    nuc_bool = 0;
end

%% basal marking
basal_bool = abs(x'-1);  % Only where x == 0;

%% Rate matrix
a = [p.Init_kmark*nuc_bool;         %nucleation
     p.kdecay*x';                   %turnover
     p.Init_basal*basal_bool];      %propogation

end