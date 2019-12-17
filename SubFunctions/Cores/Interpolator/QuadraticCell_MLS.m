function N_g = QuadraticCell_MLS(cellCount,spCount,GAUSS,CONNECT_TEMP_G,x_sp,LOCC,le,r,r_c)

% Input
% Least square matrix
Ac               = cell(cellCount,1);

for c = 1:cellCount
    Ac{c} = zeros(r_c(c),r_c(c));
end

% Function for mapping from particles to centroids using quadratic Bspline
Ng_local         = zeros(spCount,9);              % Value of quadratic bspline shape function
N_g              = cell(spCount,4);               % Value of MLS shape function

for p = 1:spCount
for c = 1:GAUSS(p)
    cpid = CONNECT_TEMP_G(p,c);
     % Compute the shape functions and gradient of the shape functions
    [Ng_local(p,c)]=Quadratic_Bspline(x_sp(p,1:2),LOCC(cpid,:),le(1),le(2));
end
 end
 
% Least square function from particle to centroids
% Compute A
for p = 1:spCount
     if r==1
         Pxy = 1;
     elseif r==3
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2)];    
     elseif r==6
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2) ; x_sp(p,1)^2 ; x_sp(p,1)*x_sp(p,2) ; x_sp(p,2)^2];   
     end
         
     PPxy = Pxy * Pxy';
     
     for c = 1:GAUSS(p)
         cpid = CONNECT_TEMP_G(p,c);
         Ac{cpid} = Ac{cpid} + Ng_local(p,c) * PPxy;
     end
end

% % Check rank D
% for c=1:cellCount
% %     if rank(A{c})<6
%     if rcond(A{c})<=9.999999e-10
%         r_c(c) = 3;
%     end
% end
 
% Compute shape function and gradients
 for p = 1:spCount
     if r==1
         Pxy = 1;
     elseif r==3
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2)];   
     elseif r==6
         Pxy = [1 ; x_sp(p,1) ; x_sp(p,2) ; x_sp(p,1)^2 ; x_sp(p,1)*x_sp(p,2) ; x_sp(p,2)^2];   
     end
     
     for c = 1:GAUSS(p)
         cpid = CONNECT_TEMP_G(p,c);
                
         if r==1
         pc = [1];
         elseif r==3
         pc = [1 ; LOCC(cpid,1) ; LOCC(cpid,2)];   
         elseif r==6
         pc = [1 ; LOCC(cpid,1) ; LOCC(cpid,2) ; LOCC(cpid,1)^2 ; LOCC(cpid,1)*LOCC(cpid,2) ; LOCC(cpid,2)^2];   
         end
         
         if r_c(cpid)==6
         
         rc = Ac{cpid} \ pc;
         N_g{p}(c) = rc' * Ng_local(p,c) * Pxy;  
         
         elseif r_c(cpid)==3
         rc = Ac{cpid}(1:3,1:3) \ pc(1:3);
         N_g{p}(c) = rc' * Ng_local(p,c) * Pxy(1:3);  
         end
         
     end
 end