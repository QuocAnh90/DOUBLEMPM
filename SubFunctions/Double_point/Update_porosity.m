function [n_sp,k_sp]=Update_porosity(n_o,F_sp,k,spCount)

n_sp                    = zeros(spCount,1);            % solid particle porousity
k_sp                    = zeros(spCount,1);            % solid particle permeability

 for spid = 1:spCount
    J = det(F_sp{spid});
    n_sp(spid) =  1 - (1 - n_o)/J;                             % porosity
    k_sp(spid) = ((1-n_o)/(1-n_sp(spid))).^2*k;
 end
