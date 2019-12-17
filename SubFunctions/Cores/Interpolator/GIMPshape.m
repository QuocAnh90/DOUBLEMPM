function [N,dN1,dN2]=GIMPshape(xp,xn,Lx,Ly,dparticle)   

lpx = dparticle(1)/2;
lpy = dparticle(2)/2;

dx = xp(1) - xn(1);
dy = xp(2) - xn(2);

if dx >= -Lx-lpx && dx <= -Lx+lpx
    Nx = (Lx+lpx+dx).^2/(4*Lx*lpx);
    dNx = (Lx+lpx+dx)/(2*Lx*lpx);
    
elseif dx > -Lx+lpx && dx <= -lpx
    Nx = 1 + dx/Lx;
    dNx = 1/Lx;
    
elseif dx > -lpx && dx <= lpx
    Nx = 1 - (dx.^2+lpx.^2)/(2*Lx*lpx);
    dNx = -dx/Lx/lpx;
    
elseif dx > lpx && dx <= Lx-lpx
    Nx = 1-dx/Lx;
    dNx = -1/Lx;
    
elseif dx > Lx-lpx && dx <= Lx+lpx
    Nx = (Lx+lpx-dx).^2/(4*Lx*lpx);
    dNx = -(Lx+lpx-dx)/(2*Lx*lpx);
    
else
    Nx = 0;
    dNx = 0;
end

if dy >= -Ly-lpy && dy <= -Ly+lpy
    Ny = (Ly+lpy+dy).^2/(4*Ly*lpy);
    dNy = (Ly+lpy+dy)/(2*Ly*lpy);
    
elseif dy > -Ly+lpy && dy <= -lpy
    Ny = 1 + dy/Ly;
    dNy = 1/Ly;
    
elseif dy > -lpy && dy <= lpy
    Ny = 1 - (dy.^2+lpy.^2)/(2*Ly*lpy);
    dNy = -dy/Ly/lpy;
    
elseif dy > lpy && dy <= Ly-lpy
    Ny = 1-dy/Ly;
    dNy = -1/Ly;
    
elseif dy > Ly-lpy && dy <= Ly+lpy
    Ny = (Ly+lpy-dy).^2/(4*Ly*lpy);
    dNy = -(Ly+lpy-dy)/(2*Ly*lpy);
    
else
    Ny = 0;
    dNy = 0;
end

    N = Nx*Ny;
    dN1 = dNx*Ny;
    dN2 = Nx*dNy;