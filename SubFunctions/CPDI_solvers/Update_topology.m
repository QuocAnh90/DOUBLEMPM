function [r1_sp,r2_sp] = Update_topology(spCount,F_sp,r1_sp,r10_sp,r2_sp,r20_sp)

for spid=1:spCount
    
    % CPDI
        r1_sp(spid,:) = (F_sp{spid} * r10_sp(spid,:)')';
        r2_sp(spid,:) = (F_sp{spid} * r20_sp(spid,:)')';    
    
%     % GIMP
%         r1_sp(spid,:) = F_sp{spid}(1,1) * r10_sp(spid,:);
%         r2_sp(spid,:) = F_sp{spid}(2,2) * r20_sp(spid,:);                   
end