function Finite_Element_wave()
    % Time and Space steps
    nt = 1000;
    nx = 250;

    % Velocity
    v = 1.;

    % Create x vectors
    dx = 1/(nx-1);
    xx = (0:(nx-1))*dx;
    x = xx;
    x(2:nx-1)=x(2:nx-1)+dx*(rand([1 nx-2])-.5)*2*0.25;

    % Number of elements
    h = diff(x);

    dt = .5*min(h)/v;

    % source function
    t = (1:nt)*dt; t0 = 3*20*dt;
    source = -1/(20*dt)^2*(t-t0).*exp(-1/(20*dt)^2*(t-t0).*(t-t0));
    source = diff(source); source(nt) = 0.;
    % source term is a spike at x=nx/2*dx; 

    f = zeros([nx 1]);f(nx/2) = 1.;
    f(1) = 0; f(nx)= 0;

    % assemble matrix Aij
    A = zeros(nx);
    for i=2:nx-1
       for j=2:nx-1
          if i==j
             A(i,j)=1/h(i-1)+1/h(i);
          elseif i==j+1
             A(i,j)=-1/h(i-1);
          elseif i+1==j
             A(i,j)=-1/h(i);
          else
             A(i,j)=0;
          end
       end
    end
     A = A*v^2;

    % assemble matrix Mij
    M=zeros(nx);
    for i=2:nx-1
       for j=2:nx-1
          if i==j
             M(i,j)=h(i-1)/3+h(i)/3;
          elseif j==i+1
             M(i,j)=h(i)/6;
          elseif j==i-1
             M(i,j)=h(i)/6;
          else
             M(i,j)=0;
          end
       end
    end

    % Solve 

    uu = zeros([nx 1]);
    u = zeros([nx nt]);


    u(2:nx-1, 1) = dt^2 * inv(M(2:nx-1,2:nx-1)) * (f(2:nx-1)*source(i)-A(2:nx-1,2:nx-1)*uu(2:nx-1)) + 2*uu(2:nx-1) - uu(2:nx-1);
    u(2:nx-1, 2) = dt^2 * inv(M(2:nx-1,2:nx-1)) * (f(2:nx-1)*source(i)-A(2:nx-1,2:nx-1)*u(2:nx-1, 1)) + 2*u(2:nx-1, 1) - uu(2:nx-1);

	% Time stepping

    for i=3:nt
	   u(2:nx-1, i) = dt^2 * inv(M(2:nx-1,2:nx-1)) * (f(2:nx-1)*source(i)-A(2:nx-1,2:nx-1)*u(2:nx-1, i-1)) + 2*u(2:nx-1, i-1) - u(2:nx-1, i-2);
    end
    
    % Create and Display 3-D graph of u(x,t)
    mesh(linspace(0, 1, nx), linspace(0, 1, nt), u');
    xlabel('x'); ylabel('t');
    title('u(x,t)');
end