function Finite_Difference_wave()
    % Length of rod
    L = 1;
    % Maxium time
    T = 1;
    
    % Number of time and space steps
    nx = 250;
    nt = 1000;
    
    % Create vectors for x and t;
    x_vec = linspace(0, L, nx);
    t_vec = linspace(0, T, nt);
    
    % Change in x and t
    dt = T/(nt - 1);
    dx = L/(nx - 1);

    % Variable
    c = 1.0;
    C = c * dt / dx; % Courant number
    
    % Boundary conditions
    bound_0 = 0.0;
    bound_nx = 0.0;
    
    % Initial condition
    %u0(1:nx) = initial_condition(nx, x_vec);%2*exp(-(x_vec-L/2).^2);
    % source function
    t = (1:nt)*dt; t0 = 3*20*dt;
    u0 = -1/(20*dt)^2*(t-t0).*exp(-1/(20*dt)^2*(t-t0).*(t-t0));
    u0 = diff(u0); u0(nt) = 0.;
    % source term is a spike at x=nx/2*dx; 
    
    % Set up u
    u = zeros ( nt, nx );
    u(1,1:nx) = u0(1:nx);
    
    % First time step
    u(2,1) = bound_0;
	u(2,2:nx-1) = C^2 * u0(3:nx)  + 2 * ( 1 - C^2 ) * u0(2:nx-1) + C^2 * u0(1:nx-2) - u0(2:nx-1);
    u(2,nx) = bound_nx;

    % Time steps loop
    for i = 3 : nt
        u(i,1) = bound_0;
        u(i,2:nx-1) = C^2 * u(i-1,3:nx)  + 2 * ( 1 - C^2 ) * u(i-1,2:nx-1) + C^2 * u(i-1,1:nx-2) - u(i-2,2:nx-1);
        u(i,nx) = bound_nx;
    end
	
	% Display 3-D graph of u(x,t)
% 	[X, T] = meshgrid(x_vec, t_vec);
% 	surfl(X', T', u');
% 	xlabel('x'); ylabel('t');
% 	title('Solution u(x,t)');
    mesh(x_vec, t_vec, u);
    xlabel('x'); ylabel('t');
    title('u(x,t)');
end