% Dimensional Splitting and High Resolution Method
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
   
   q1 = zeros(N,N); q2 = zeros(N,N); q3 = zeros(N,N);
   theta = zeros(N,N); limiter = zeros(N,N); flux = zeros(2*N+1,2*N+1); fluxp = zeros(2*N+1,2*N+1);
   flm = zeros(2*N+1,2*N+1);  fhm = zeros(2*N+1,2*N+1); 
   q = q0;
   
   figure;
   % Create the video write with 1 fps
    writerObj = VideoWriter('strang_minmod.avi');
    % Set the seconds per image
    writerObj.FrameRate = 20;
    % Open the video writer
    open(writerObj);
   
for t = 1:1:T   % Loop over time      
    %Looping over vertical edges (i-1/2,j)
    for j = 2:2:2*N+1
        for i = 1:2:2*N-1  
             %compute the flux function          
             if i == 1
               flm(i,j) = uhalfm(i,j)*q(1,j/2)+uhalfp(i,j)*q(N,j/2) ;   %upwinding flux
               fhm(i,j) = uhalfm(i,j)*q(1,j/2)+uhalfp(i,j)*q(N,j/2)...  %lax-wendroff flux
                    +(0.5*abs(uhalf(i,j))*(1-(dt/dx)*abs(uhalf(i,j)))*(q((i+1)/2,j/2)-q(N,j/2)));
             else 
               flm(i,j) = uhalfm(i,j)*q((i+1)/2,j/2)+uhalfp(i,j)*q((i-1)/2,j/2) ;   %upwinding flux
               fhm(i,j) = uhalfm(i,j)*q((i+1)/2,j/2)+uhalfp(i,j)*q((i-1)/2,j/2)...  %lax-wendroff flux
                    +(0.5*abs(uhalf(i,j))*(1-(dt/dx)*abs(uhalf(i,j)))*(q((i+1)/2,j/2)-q((i-1)/2,j/2))); 
             end
      
             %smoothness indicator 
             if  i == 1 || i == 2*N-1
                 theta((i+1)/2,j/2) = (q((i+1)/2,j/2)-q(N,j/2))/(q(1,j/2)-q((i+1)/2,j/2));
             else
                 theta((i+1)/2,j/2) = (q((i+1)/2,j/2)-q((i-1)/2,j/2))/(q((i+3)/2,j/2)-q((i+1)/2,j/2));
             end 
    
             % Minmod limiter
            limiter((i+1)/2,j/2) = min(1.0,theta((i+1)/2,j/2));
            limiter((i+1)/2,j/2) = max(0,limiter((i+1)/2,j/2));
             
             %flux at the left edge (i-1/2)
             flux(i,j) = flm(i,j)+limiter((i+1)/2,j/2)*(fhm(i,j)-flm(i,j));
        end 
    end
    
     for i = 1:1:N
         for j = 1:1:N   
             q1(i,j) = q(i,j)-(0.5*dt/dx)*(flux(2*i+1,2*j)-flux(2*i-1,2*j));
         end
     end
     
    % Looping over horizontal edges (i,j-1/2)
    for j = 1:2:2*N-1   % to avoid duplicates in the last cell
        for i = 2:2:2*N+1        

             %compute the low res flux function     (FL(i,j-1/2))      
             if j == 1 
               flm(i,j) = vhalfm(i,j)*q1(i/2,1)+vhalfp(i,j)*q1(i/2,N) ;   %upwinding flux
               fhm(i,j) = vhalfm(i,j)*q1(i/2,1)+vhalfp(i,j)*q1(i/2,N)...  %lax-wendroff flux
                    +(0.5*abs(vhalf(i,j))*(1-(dt/dx)*abs(vhalf(i,j)))*(q1(i/2,(j+1)/2)-q1(i/2,N)));
             else 
               flm(i,j) = vhalfm(i,j)*q1(i/2,(j+1)/2)+vhalfp(i,j)*q1(i/2,(j-1)/2);   %upwinding flux
               fhm(i,j) = vhalfm(i,j)*q1(i/2,(j+1)/2)+vhalfp(i,j)*q1(i/2,(j-1)/2)...  %lax-wendroff flux
                    +(0.5*abs(vhalf(i,j))*(1-(dt/dx)*abs(vhalf(i,j)))*(q1(i/2,(j+1)/2)-q1(i/2,(j-1)/2))); 
             end
      
             %smoothness indicator 
             if  j == 1 || j == 2*N-1
                 theta(i/2,(j+1)/2) = (q1(i/2,(j+1)/2)-q1(i/2,N))/(q1(i/2,1)-q1(i/2,(j+1)/2));
             else
                 theta(i/2,(j+1)/2) = (q1(i/2,(j+1)/2)-q1(i/2,(j-1)/2))/(q1(i/2,(j+3)/2)-q1(i/2,(j+1)/2));
             end 
    
             % Minmod limiter
            limiter(i/2,(j+1)/2) = min(1.0,theta(i/2,(j+1)/2));
            limiter(i/2,(j+1)/2) = max(0,limiter(i/2,(j+1)/2));
             
             %flux at the bottom edge (j-1/2)
             flux(i,j) = flm(i,j)+limiter(i/2,(j+1)/2)*(fhm(i,j)-flm(i,j));
        end
    end
    
     for i = 1:1:N
         for j = 1:1:N   
             q2(i,j) = q1(i,j)-(dt/dy)*(flux(2*i,2*j+1)-flux(2*i,2*j-1));
         end
     end 
     
    % Looping over vertical edges (i-1/2,j)
    for j = 2:2:2*N+1
        for i = 1:2:2*N-1  
             %compute the flux function          
             if i == 1
               flm(i,j) = uhalfm(i,j)*q2(1,j/2)+uhalfp(i,j)*q2(N,j/2) ;   %upwinding flux
               fhm(i,j) = uhalfm(i,j)*q2(1,j/2)+uhalfp(i,j)*q2(N,j/2)...  %lax-wendroff flux
                    +(0.5*abs(uhalf(i,j))*(1-(dt/dx)*abs(uhalf(i,j)))*(q2((i+1)/2,j/2)-q2(N,j/2)));
             else 
               flm(i,j) = uhalfm(i,j)*q2((i+1)/2,j/2)+uhalfp(i,j)*q2((i-1)/2,j/2) ;   %upwinding flux
               fhm(i,j) = uhalfm(i,j)*q2((i+1)/2,j/2)+uhalfp(i,j)*q2((i-1)/2,j/2)...  %lax-wendroff flux
                    +(0.5*abs(uhalf(i,j))*(1-(dt/dx)*abs(uhalf(i,j)))*(q2((i+1)/2,j/2)-q2((i-1)/2,j/2))); 
             end
      
             %smoothness indicator 
             if  i == 1 || i == 2*N-1
                 theta((i+1)/2,j/2) = (q2((i+1)/2,j/2)-q2(N,j/2))/(q2(1,j/2)-q2((i+1)/2,j/2));
             else
                 theta((i+1)/2,j/2) = (q2((i+1)/2,j/2)-q2((i-1)/2,j/2))/(q2((i+3)/2,j/2)-q2((i+1)/2,j/2));
             end 
    
             % Minmod limiter
            limiter((i+1)/2,j/2) = min(1.0,theta((i+1)/2,j/2));
            limiter((i+1)/2,j/2) = max(0,limiter((i+1)/2,j/2));
             
             %flux at the left edge (i-1/2)
             flux(i,j) = flm(i,j)+limiter((i+1)/2,j/2)*(fhm(i,j)-flm(i,j));
        end 
    end
    
     for i = 1:1:N
         for j = 1:1:N   
             q3(i,j) = q2(i,j)-(0.5*dt/dx)*(flux(2*i+1,2*j)-flux(2*i-1,2*j));
         end
     end
     
      q = q3; 
      surf(xp,yp,q3,'LineStyle',':');
      colormap(jet(1000));
      colorbar;
      title(sprintf('Dimensional Splitting and High Resolution solution Scheme (N=%d,t=%3.4f)', t,t*dt));
      zlim([-0.2 1.2]); 
      xlabel('x');ylabel('y'); zlabel('q');
      drawnow; 
      F(t) = getframe(gcf);
      frame = F(t);
      writeVideo(writerObj,frame);
      
end
    
    % Plot of 2D contour plot
    figure;
    % Define a vector to show certain, important contour lines only
    v = [0.08 0.18 0.32 0.44 0.56 0.65 0.77 0.82 0.94 0.102 0.2174 0.3358 0.4284 0.6059 0.69];
    contour(xp,yp,q3,v);
    axis equal;
    xlabel('x');ylabel('y'); 