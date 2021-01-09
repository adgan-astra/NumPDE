% Corner Transport Upwinding Scheme
clear workspace;

% Define the bounds of the domain
ax = -1; ay = -1;
bx = 1; by = 1;

% grid spacing
dx = 1/64;
dy = 1/64;

% Step size
dt = 0.4*dx;   % CTU is unstable if ratio is greater than 0.5
T = ceil(pi/dt);

% Number of grid points is the same in x and y
N = (bx-ax)/dx;

% Define a grid with cell centers and edge centers
x = zeros(1,N); y = zeros(1,N);
for i = 1:1:(2*N+1)
    x(i) = ax + 0.5*(i-1)*dx;
    y(i) = ay +0.5*(i-1)*dy;
end

% Initial Conditions
q0 = zeros(N,N); 

for i = 1:1:N
    for j = 1:1:N
        if sqrt((x(2*i)+0.45)^2+y(2*j)^2) < 0.35
             q0(i,j) = 1-(sqrt((x(2*i)+0.45)^2+y(2*j)^2))/0.35;
        end
    end
end

for i = 1:1:N
    for j = 1:1:N
        if x(2*i) < 0.6 && x(2*i) > 0.1 && y(2*j) > -0.25  && y(2*j) < 0.25
             q0(i,j) = 1;
        end
    end
end
q = q0;

% Initial/ analytical solution
% surf(xp,yp,q0);xlim([-1,1]); ylim([-1 1]);
% zlim([0 1]); 
% xlabel('x');ylabel('y'); zlabel('q');

% coordinates of cell centers
xp = zeros(1,N); yp = zeros(1,N);
for i = 1:1:N
    xp(i) = x(2*i);
    yp(i) = y(2*i);
end

% Define edge centered velocities u and v
uhalf = zeros(2*N+1, 2*N+1); vhalf = zeros(2*N+1,2*N+1);
uhalfp = zeros(2*N+1, 2*N+1); uhalfm = zeros(2*N+1, 2*N+1);
vhalfp = zeros(2*N+1, 2*N+1); vhalfm = zeros(2*N+1, 2*N+1);

 for i = 1:2:2*N+1
    for j = 2:2:2*N+1
        % Now take the velocity averages at each vertical edge
        %u(i-1/2,j)
        uhalf(i,j) = 2*y(j);
        uhalfp(i,j) = max(0,uhalf(i,j));
        uhalfm(i,j) = min(0,uhalf(i,j));
     end
end

   for j = 1:2:2*N+1   
       for i = 2:2:2*N+1
            % Now take the velocity averages at each horizontal edge
            %v(i,j-1/2)   
            vhalf(i,j) = -2*x(i);
            vhalfp(i,j) = max(0,vhalf(i,j));
            vhalfm(i,j) = min(0,vhalf(i,j));
       end
   end

% Corrective G fluxes
G1 = zeros(2*N+1,2*N+1); G2 = zeros(2*N+1,2*N+1); 
G3 = zeros(2*N+1,2*N+1); G4 = zeros(2*N+1,2*N+1);
G = zeros(2*N+1,2*N+1);   

% Corrective F fluxes
F1 = zeros(2*N+1,2*N+1); F2 = zeros(2*N+1,2*N+1); 
F3 = zeros(2*N+1,2*N+1); F4 = zeros(2*N+1,2*N+1);
F = zeros(2*N+1,2*N+1);

qnew = zeros(N,N);  
figure;

for t = 1:1:T   % Loop over time
        
    for i = 1:1:N
      for j = 1:1:N
         if i == 1 || i == N || j == 1 || j == N
             % Use periodic boundary conditions
           qnew(i,j) = q(i,j) - (dt/dx)*(uhalfp(2*i-1,2*j)*(q(i,j)-q(N,j))+uhalfm(2*i+1,2*j)*(q(1,j)-q(i,j)))...
                  -(dt/dy)*(vhalfp(2*i,2*j-1)*(q(i,j)-q(i,N))+vhalfm(2*i,2*j+1)*(q(i,1)-q(i,j)))...
                  -(dt/dy)*(G(2*i,2*N+1)-G(2*i,2*N-1))...
                -(dt/dx)*(F(2*N+1,2*j)-F(2*N-1,2*j));
         else
           qnew(i,j) = q(i,j) - (dt/dx)*(uhalfp(2*i-1,2*j)*(q(i,j)-q(i-1,j))+uhalfm(2*i+1,2*j)*(q(i+1,j)-q(i,j)))...
                  -(dt/dy)*(vhalfp(2*i,2*j-1)*(q(i,j)-q(i,j-1))+vhalfm(2*i,2*j+1)*(q(i,j+1)-q(i,j)))...
                  -(dt/dy)*(G(2*i,2*j+1)-G(2*i,2*j-1))...
                  -(dt/dx)*(F(2*i+1,2*j)-F(2*i-1,2*j)); 
         end
      end
    end
    
    q = qnew;     
  
    % Looping over vertical edges (i-1/2,j)
    for j = 2:2:2*N+1
        for i = 1:2:2*N-1  % Done to avoid duplicates in the last cell
                if i == 1
                    G1(2*N,j-1) = -0.5*uhalfm(i,j)*vhalfm(2*N,j-1)*dt/dx*(q(1,j/2)-q(N,j/2));
                    G2(2*N,j+1) = -0.5*uhalfm(i,j)*vhalfp(2*N,j+1)*dt/dx*(q(1,j/2)-q(N,j/2));
                    G3(2,j-1) = -0.5*uhalfp(i,j)*vhalfm(2,j-1)*dt/dx*(q(1,j/2)-q(N,j/2));
                    G4(2,j+1) = -0.5*uhalfp(i,j)*vhalfp(2,j+1)*dt/dx*(q(1,j/2)-q(N,j/2));    
                else   
                    G1(i-1,j-1) =  -0.5*uhalfm(i,j)*vhalfm(i-1,j-1)*dt/dx*(q((i+1)/2,j/2)-q((i-1)/2,j/2)); % Last term is at j= 2N-1
                    G2(i-1,j+1) =  -0.5*uhalfm(i,j)*vhalfp(i-1,j+1)*dt/dx*(q((i+1)/2,j/2)-q((i-1)/2,j/2));  
                    G3(i+1,j-1) =  -0.5*uhalfp(i,j)*vhalfm(i+1,j-1)*dt/dx*(q((i+1)/2,j/2)-q((i-1)/2,j/2)); % Last term is at j= 2N-1
                    G4(i+1,j+1) =  -0.5*uhalfp(i,j)*vhalfp(i+1,j+1)*dt/dx*(q((i+1)/2,j/2)-q((i-1)/2,j/2));
                end
        end
    end
  
     % Accumulate fluxes for every horizontal interface by summing G1,G2,G3,G4
     % appropriately
     for j = 2:2:2*N
         for i = 1:2:2*N-1
              if j == 2*N 
                  G(i+1,2*N+1) =G1(i+1,1)+G2(i+1,2*N+1)+G3(i+1,1)+G4(i+1,2*N+1);   % Upper boundary
              else 
                  G(i+1,j+1) = G1(i+1,j+1)+G2(i+1,j+1)+G3(i+1,j+1)+G4(i+1,j+1);
              end
         end
     end

    % Looping over horizontal edges (i,j-1/2)
       for j = 1:2:2*N-1   % to avoid duplicates in the last cell
           for i = 2:2:2*N+1
                if j == 1 
                    F1(i-1,2*N) = -0.5*vhalfm(i,j)*uhalfm(i-1,2*N)*dt/dy*(q(i/2,1)-q(i/2,N));
                    F2(i-1,2) =  -0.5*vhalfp(i,j)*uhalfm(i-1,2)*dt/dy*(q(i/2,1)-q(i/2,N));
                    F3(i+1,2*N) = -0.5*vhalfm(i,j)*uhalfp(i+1,2*N)*dt/dy*(q(i/2,1)-q(i/2,N));
                    F4(i+1,2) = -0.5*vhalfp(i,j)*uhalfp(i+1,2)*dt/dy*(q(i/2,1)-q(i/2,N));    
                else   
                    F1(i-1,j-1) = -0.5*vhalfm(i,j)*uhalfm(i-1,j-1)*dt/dy*(q(i/2,(j+1)/2)-q(i/2,(j-1)/2));  % Last term is at i = 2N-1
                    F2(i-1,j+1) = -0.5*vhalfp(i,j)*uhalfm(i-1,j+1)*dt/dy*(q(i/2,(j+1)/2)-q(i/2,(j-1)/2));  % Last term is at i = 2N-1
                    F3(i+1,j-1) = -0.5*vhalfm(i,j)*uhalfp(i+1,j-1)*dt/dy*(q(i/2,(j+1)/2)-q(i/2,(j-1)/2));   
                    F4(i+1,j+1) = -0.5*vhalfp(i,j)*uhalfp(i+1,j+1)*dt/dy*(q(i/2,(j+1)/2)-q(i/2,(j-1)/2));
                end
           end
       end

    % Accumulate fluxes for every vertical edge by summing F1,F2,F3,F4
    % appropriately
    for j = 1:2:2*N-1
       for i = 2:2:2*N
              if i == 2*N
                  F(2*N+1,j+1) = F3(2*N+1,j+1)+F4(2*N+1,j+1)+F1(1,j+1)+F2(1,j+1); % Right boundary
              else 
                  F(i+1,j+1) = F1(i+1,j+1)+F2(i+1,j+1)+F3(i+1,j+1)+F4(i+1,j+1);
              end
       end    
    end
     
  surf(xp,yp,qnew,'LineStyle',':');
  colormap(jet(1000));
  colorbar;
  title(sprintf('CTU solution (N=%d,t=%3.4f)', t,t*dt));
  zlim([0 1]); 
  xlabel('x');ylabel('y'); zlabel('q');
  drawnow; pause(0.000001);     
end

% Plot of 2D contour plot
figure;
% Define a vector to show certain, important contour lines only
v = [0.08 0.28 0.35 0.47 0.58 0.99 0.87 0.21 0.32 0.44 0.55 0.61 0.69];
contour(xp,yp,qnew,v);
xlabel('x'); ylabel('y');
axis equal;

% Plot of difference in analytical and CTU solution
figure;
g = q0-qnew;
surf(xp,yp,g);
colormap(jet(1000));
colorbar;
xlabel('x'); ylabel('y'); zlabel('q_{analytical}-q_{CTU}');