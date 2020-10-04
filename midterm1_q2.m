% Midterm 1 - Q2
clear all;
close all;

dx = 0.01;
dt = 4*dx;
a = 0.0;
B = 1.0;
N = (B-a)/dx;   % Spatial nodes
x = linspace(a+0.5*dx,B-0.5*dx,N);  % Include the ghost cells
u = zeros(1,N);

% Define start and endtimes
StartTime = 0.0;
EndTime = 4.0;

% compute number of time steps
Nsteps = ceil((EndTime-StartTime)/dt);

% initial condition
for i=1:1:N
 u(i) = exactsolution(x(i));
end

t = StartTime;

% Initialize coefficient matrix A
A = zeros(N,N);

% Start populating A:

A(1,N) = -1/dt;     % Due to boundary condition
A(N,1) = 1/dt;      % Due to boundary condition
A(N,N) = 1/dt;

% Populate the rest of matrix A
for i = 1:N-1
        A(i+1,i) = -1/dt;
        A(i,i) = 1/dt;
        A(i,i+1) = 1/dt;
end 
    
b = zeros(1,N);
% Populate RHS vector b:
b(1) = (1/dt)*u(N)+(1/dt)*u(1)-(1/dt)*u(2);
b(N) = (1/dt)*u(N-1)+(1/dt)*u(N)-(1/dt)*u(1);
for i = 2:N-1
    b(i) = (1/dt)*u(i-1)+(1/dt)*u(i)-(1/dt)*u(i+1);
end
b = transpose(b);

% Loop through time
for k = 1:Nsteps
    
    unew = A\b;  % Update solution after each time-step
    u = unew;
    t  = t+dt;   % update time
            
    % Re-populate RHS vector b:
    b(1) = (1/dt)*u(N)+(1/dt)*u(1)-(1/dt)*u(2);
    b(N) = (1/dt)*u(N-1)+(1/dt)*u(N)-(1/dt)*u(1);
    for i = 2:N-1
        b(i) = (1/dt)*u(i-1)+(1/dt)*u(i)-(1/dt)*u(i+1);
    end
        
     plot(x, u, 'bo-', 'LineWidth', 1.2); % plot numerical solution
     xlabel('x','fontsize',14);
     ylabel('u','fontsize',14);
     hold on;
     xplot = linspace(a, B, 1000); % plot exact solution at lots of points
     xexact = xplot-t ;            % coordinate for exact pulse
     uexact = zeros(1,1000);
     for i=1:1:1000
       uexact(i) = exactsolution(xexact(i));
     end
     plot(xplot, uexact, 'r-', 'LineWidth', 2);
     hold off;
     axis([a B -.2 1.2]);
     legend(sprintf('Numerical solution (N=%d,t=%3.4f)', k,t), 'Exact solution');
     drawnow; pause(0.02);
     
end   