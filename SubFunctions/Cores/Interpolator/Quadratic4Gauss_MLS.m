function N_g = Quadratic4Gauss_MLS(cellCount,spCount,GAUSS,CONNECT_TEMP_G,x_sp,GaussCount,x_gauss,y_gauss,le,r,r_c)

% Input
% Least square matrix
Ac               = cell(cellCount,GaussCount);

% Function for mapping from particles to centroids using quadratic Bspline
Ng_local         = cell(spCount,1);                     % Value of quadratic bspline shape function
N_g              = cell(spCount,1);                     % Value of MLS shape function

for c = 1:cellCount
for g = 1:GaussCount
    Ac{c,g} = zeros(r_c(c),r_c(c));
end
end

for p = 1:spCount
    Ng_local{p} = zeros(GAUSS(p),GaussCount);
    N_g{p} = zeros(GAUSS(p),GaussCount);
end

% Compute the quadratic Bspline shape function for Gauss point
for p = 1:spCount
for c = 1:GAUSS(p)
    cpid = CONNECT_TEMP_G(p,c);

for g = 1:GaussCount
     % Compute the shape functions and gradient of the shape functions
        position_Gauss = [x_gauss(cpid,g);y_gauss(cpid,g)];
        Ng_local{p}(c,g)=Quadratic_Bspline(x_sp(p,1:2),position_Gauss,le(1),le(2));
%         Ng_local{p}(c,g)=Cubic_Bspline(x_sp(p,1:2),position_Gauss,le(1),le(2));
 end
end
end

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
     for g = 1:GaussCount

         Ac{cpid,g} = Ac{cpid,g} + Ng_local{p}(c,g) * PPxy;
     end
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
         
         for g = 1:GaussCount
         position_Gauss = [x_gauss(cpid,g);y_gauss(cpid,g)];

         if r==1
         pc = [1];
         elseif r==3
         pc = [1 ; position_Gauss(1) ; position_Gauss(2)];   
         elseif r==6
         pc = [1 ; position_Gauss(1) ; position_Gauss(2) ; position_Gauss(1)^2 ; position_Gauss(1)*position_Gauss(2) ; position_Gauss(2)^2];   
         end
         
         if r_c(cpid)==6
         
         rc = Ac{cpid,g} \ pc;
         N_g{p}(c,g) = rc' * Ng_local{p}(c,g) * Pxy;  
         
         elseif r_c(cpid)==3
         rc = Ac{cpid,g}(1:3,1:3) \ pc(1:3);
         N_g{p}(c,g) = rc' * Ng_local{p}(c,g) * Pxy(1:3);  
         end
         end
     end
 end