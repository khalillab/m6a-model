function a = lattice_clus_withRW(x, p)
% Return reaction propensities given current state x
% x = current state of lattice
% p = structure array with rate parameters

L = length(x);

%% nucleation
nuc_bool = abs(x(1)-1); % When x(1) == 0;

%% basal marking
basal_bool = abs(x'-1);  % Only where x == 0;

%% propogation logic
prop_bool = zeros(L,1);

    % first site (edge case)
    if (x(2)==1) && (x(1)==0)
        prop_bool(1) = 1;
    end

    % last site (edge case)
    if (x(end-1)==1) && (x(end)==0)
        prop_bool(end) = 1;
    end

    % all other sites
%     for i=2:(L-1)
%         if (x(i)==0) && ((x(i-1)==1)||(x(i+1)==1))
%             prop_bool(i) = 1;
%         end
%     end

    % all other sites (modified for speed)
    b = basal_bool(2:L-1);
    prop_bool(2:L-1) = (b>basal_bool(1:L-2))|(b>basal_bool(3:L));

    
%% rate matrix
a = cat(1, ...
        p.Init_kmark*nuc_bool, ...                  % nucleation
        (p.Init_prop+p.RW_prop)*prop_bool, ...      % propogation
        p.kdecay*x', ...                            % turnover
        (p.Init_basal+p.RW_basal)*basal_bool);      % basal marking
    
end