function  Assignment3_Part2
clear
width = 1;
len = 2;
V0 = .1;
dx = 0.025; % Mesh spacing along x
dy = 0.025; % Mesh spacing along y
nx = ceil(len/dx); % Number of points along x
ny = ceil(width/dy); % Number of points along y
dx = len/nx;
dy = width/ny;
Lb = 0.4; % Length of the highly resistive regions
Wb = 0.4; % Width of the highly resistive regions
sigma = 1; % Nominal conductivity
sigma_insulating = 1e-2; % Conductivity in highly resistive regions

% Construct the C matrix:
C = sigma.*ones(ny,nx);
sigma_map = zeros(ny,nx);

for x=1:nx
    for y=1:ny
        xx = x*dx;
        yy = y*dy;
        
        % The resistivity is made high in the rectangular regions:
        if(xx <= (len+Lb)/2 && xx >= (len-Lb)/2 && (yy >= width-Wb || yy <= Wb))
            sigma_map(y,x) = sigma-sigma_insulating;
        end
    end
end


sigma_map = imgaussfilt(sigma_map, 1);
C = C - sigma_map;

%%


figure(5);
surf(C,'EdgeColor','none');
colorbar
view(0,90)
title('Conductivity');
xlim([0 ceil(len/dx)])
ylim([0 ceil(width/dy)])
xlabel('x (m)');
ylabel('y (m)');
grid on;

%%
G = zeros(nx*ny,nx*ny);
F = zeros(nx*ny,1);

dx2 = 1./(dx.^2);
dy2 = 1./(dy.^2);

for x=2:(nx-1)
    for y=2:(ny-1)
        index = mapCoordinate(x,y,nx);
        
        % Apply the equation derived earlier:
        G(index,index) = -2.*C(y,x).*(dx2 + dy2);
        G(index, mapCoordinate(x+1,y,nx)) = dx2.*(0.25.*(C(y,x+1) - C(y,x-1)) + C(y,x));
        G(index, mapCoordinate(x-1,y,nx)) = dx2.*(-0.25.*(C(y,x+1) - C(y,x-1)) + C(y,x));
        
        G(index, mapCoordinate(x,y+1,nx)) = dy2.*(0.25.*(C(y+1,x) - C(y-1,x)) + C(y,x));
        G(index, mapCoordinate(x,y-1,nx)) = dy2.*(-0.25.*(C(y+1,x) - C(y-1,x)) + C(y,x));
    end
end
for x=2:(nx-1)
    index = mapCoordinate(x,1,nx);
    G(index,index) = 1;
    G(index,mapCoordinate(x,2,nx)) = -1;
    F(index) = 0;
    
    index = mapCoordinate(x,ny,nx);
    G(index,index) = 1;
    G(index,mapCoordinate(x,ny-1,nx)) = -1;
    F(index) = 0;
end

% The vertical boundaries
for y=1:ny
    index = mapCoordinate(1,y,nx);
    G(index,index) = 1;
    F(index) = V0;
    
    index = mapCoordinate(nx,y,nx);
    G(index,index) = 1;
    F(index) = 0;
end

V = G\F;
V = reshape(V,[],ny)';

%%

figure(6);
contourf(V,30);
colorbar
view(0,90);
xlabel('x (m)');
ylabel('y (m)');
xlim([0 ceil(len/dx)])
ylim([0 ceil(width/dy)])
title('Electric Potential (V)');

figure(7);
[Ex,Ey] = gradient(-V);
quiver(Ex,Ey,5);
xlabel('x (m)');
ylabel('y (m)');
title('Electric Field (V/m)');
xlim([0 ceil(len/dx)])
ylim([0 ceil(width/dy)])
grid on;
end
function index = mapCoordinate(x,y,xn)

index = (y-1).*xn + x;

end

