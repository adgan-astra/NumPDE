% Midterm 1 - Q1

clear workspace;
close all;

% Define variables
a = 0.0;
b = 1.0;
dx = 0.01;
N = (b-a)/dx;   % Spatial nodes
dt = 0.75*dx;   % time step
v = 1;          % velocity = 1

x = linspace(a+0.5*dx,b-0.5*dx,N);  % Include the ghost cells
u = zeros(1,N);

% start time
StartTime = 0.0;
EndTime = 2.0;

% compute number of time steps
Nsteps = ceil((EndTime-StartTime)/dt);

% initial condition
for i=1:1:N
     u(i) = exactsolution(x(i)); % First I.C
   % u(i) = exactsolution1(x(i)); % Second I.C

end

t = StartTime;

% Loop through time
for n = 1:1:Nsteps;
     if n == (Nsteps);   % To make sure loop stops at time = 2
         dt = 0.005;
     end
     
     % Compute the flux function
     f = v*u;    
     
     %flux at the right edge (j+1/2) 
     flux = f;
     
     %periodic condition: flux at the left edge (j-1/2) 
     fluxm1 = [flux(N),flux(1:N-1)];
     
     % Evaluate u_n+1
     unew = u-v*(dt/dx)*(flux-fluxm1);
     u = unew;    % finish one time step
     t  = t+dt;   % update time
     
     % Plot the exact and numerical solutions
     plot(x, u, 'bo-', 'LineWidth', 1.2); % plot numerical solution
     hold on;
     xplot = linspace(a, b, 1000); % plot exact solution at lots of points
     xexact = xplot-t ;            % coordinate for exact pulse
     uexact = zeros(1,1000);
     for i=1:1:1000
       uexact(i) = exactsolution(xexact(i));
       %uexact(i) = exactsolution1(xexact(i));  
     end
     plot(xplot, uexact, 'r-', 'LineWidth', 2);
     xlabel('x','fontsize',14);
     ylabel('u','fontsize',14);
     hold off;
     axis([a b -.2 1.2]);
     legend(sprintf('Numerical solution (N=%d,t=%3.4f)', n,t), 'Exact solution');
     drawnow; pause(0.02);
    
end
