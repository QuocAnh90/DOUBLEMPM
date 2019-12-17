function [n_wp,k_wp,V_wp]=Update_liquid_porosity_particle(NODES,CONNECT,N,wpCount,k_si,n_si,nmass_si,n_o,V_wpo)

n_wp                    = zeros(wpCount,1);                 % liquid particle porosity
k_wp                    = zeros(wpCount,1);                 % liquid particle permeability

for wp = 1:wpCount
     for j = 1:NODES(wp)
         npid                           = CONNECT{wp}(j);
              if nmass_si(npid)==0
                continue
              end         
         n_wp(wp)           = n_wp(wp) + n_si(npid) * N{wp}(j);
         k_wp(wp)           = k_wp(wp) + k_si(npid) * N{wp}(j);
     end   
end

for wp=1:wpCount
    V_wp(wp) = (n_o * V_wpo(wp))/n_wp(wp);
end
