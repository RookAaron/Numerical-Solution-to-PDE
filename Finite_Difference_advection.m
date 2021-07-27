function Finite_Difference_advection()
    % Length of rod
    L = 1;
    % Maximum time
    T = 1;
    
    % Number of time and space steps
    nx = 151;
    nt = 1000;
    
    % Create vectors for x and t;
    x_vec = linspace(0, L, nx);
    t_vec = linspace(0, T, nt);
    
    % Change in x and t
    dt = T/nt;
    dx = L/(nx - 1);

    % Variable
    c = 1.0;
    C = c * dt / dx % Courant number
    
    % Initial condition
    u0(1:nx) = initial_condition ( nx, x_vec );
    
    % Set up u
    u = zeros ( nt, nx );
    u(1,1:nx) = u0(1:nx);
        
    % Space order for Euler Equation
    backward_space = [nx, 1:nx-2, nx-1];
    space = [1, 2:nx-1, nx];
    forward_space = [2, 3:nx, 1];
    
    % Loop to calculate all values
    for j = 2:nt
        u(j,:) = u(j-1,space) - C/2*(u(j-1,forward_space)-u(j-1,backward_space));
    end
	
	% Create and Display 3-D graph of u(x,t)
    mesh(x_vec, t_vec, u);
    xlabel('x'); ylabel('t');
    title('u(x,t)');
end

function u0 = initial_condition(nx, x_vec)
    u0 = zeros(1,nx);
    for i=1:nx
        if 0.4 <= x_vec(i) && x_vec(i) <= 0.6
            u0(i) = (10.0*x_vec(i)-4)^2 *(6-10*x_vec(i))^2;
        end
    end
end