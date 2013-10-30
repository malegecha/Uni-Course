

function A2matlab()

    close all; clear all; clc;

    global T M K s rho Qcpu C k Tair h N_p N_f N_e N_b ...
           Points Faces Elements Boundaries Delta_t l;

    [Points, Faces, Elements, Boundaries, N_p, N_f, N_e, N_b] = Box();
    
    % Simulation parameters
    x_min           =  0.00;
    x_max           =  0.02;
    y_min           =  0.00;
    y_max           =  0.02;
    t_min           =  0.00;
    t_max           =  100.00;        
    Delta_t         =  1;
    t               =  t_min:Delta_t:t_max;
    N_t             =  length(t);
    rho             =  8954;
    C               =  380;
    k               =  386;
    Tair            =  300;
    h               =  100;
    Qcpu            =  40000;
	
    % Allocate arrays
    T               = zeros (N_p, N_t);
    M               = sparse(N_p, N_p);
    K               = sparse(N_p, N_p);
    s               = sparse(N_p, 1);
    % Set initial condition
    T(:,1)          = Tair;

    assembleSystem();
    % Implicit Euler method for the spatial discretization
    A               = M - Delta_t*K;
    Free            = 1:N_p;

    % Set up animated visualization of results
    figure('WindowStyle', 'docked');
    Solution  = trisurf(Faces(Boundaries(1).indices,:), Points(:,1),...
                Points(:,2), Points(:,3),T(:,1)); hold on;
    Solution2 = trisurf(Faces(Boundaries(2).indices,:), Points(:,1),...
                Points(:,2), Points(:,3),T(:,1)); hold on;
    axis( [x_min x_max y_min y_max]);
    grid on;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    view([45 25]);
    colormap('Winter');
    drawnow;
         
    % Main time marching loop
    for l=1:N_t-1
        b        	= M*T(:,l) + Delta_t*s;
        T(Free,l+1)	= A(Free,Free)\b(Free);
        % Plot T at every timestep
          set(Solution,'CData', T(:,l+1));
          set(Solution2, 'CData', T(:,l+1));
        title(['t = ' num2str(t(l+1))]);
        drawnow;
        
    end
    figure;
    plot(t,T(1,:));
    title('Change of Temperature')
    axis([t_min t_max 300 345]);
    xlabel('time');
    ylabel('temperature');
    grid on;
end
    
function assembleSystem()

    global M K k s h C rho Tair Qcpu N_f N_e N_b Points Faces ...
           Elements Boundaries Gamma Omega ;
    % assemble system
	Omega       = zeros(N_e, 1);
    Gamma       = zeros(N_f, 1);
	M_e         = [2, 1, 1, 1; 
                   1, 2, 1, 1; 
                   1, 1, 2, 1;
                   1, 1, 1, 2];
    K_f         =[2, 1, 1;
                  1, 2, 1;
                  1, 1, 2];
 	s_e         = [1;
                   1;
                   1];
    % Calculate face lengths
    for f=1:N_f
        x           = Points(Faces(f, :), 1);
        y           = Points(Faces(f, :), 2);
        z           = Points(Faces(f, :), 3);
        Gamma(f)	= sqrt(((y(2)-y(1))*(z(3)-z(1)) ... 
                      - (z(2)-z(1))*(y(3)-y(1)))^2 ...
                      + ((z(2)-z(1))*(x(3)-x(1)) ...
                      - (x(2)-x(1))*(z(3)-z(1)))^2 ...
                      + ((x(2)-x(1))*(y(3)-y(1)) ...
                      - (y(2)-y(1))*(x(3)-x(1)))^2)/2;
    end
    % Calculate element areas
    for e=1:N_e	
        x           = Points(Elements(e, :), 1);
        y           = Points(Elements(e, :), 2);
        z           = Points(Elements(e, :), 3);
        Omega(e)	= abs( x(1)*y(2)*z(3) - x(1)*y(3)*z(2) ...
                         - x(2)*y(1)*z(3) + x(2)*y(3)*z(1) ...
                         + x(3)*y(1)*z(2) - x(3)*y(2)*z(1) ...
                         - x(1)*y(2)*z(4) + x(1)*y(4)*z(2) ...
                         + x(2)*y(1)*z(4) - x(2)*y(4)*z(1) ...
                         - x(4)*y(1)*z(2) + x(4)*y(2)*z(1) ...
                         + x(1)*y(3)*z(4) - x(1)*y(4)*z(3) ...
                         - x(3)*y(1)*z(4) + x(3)*y(4)*z(1) ...
                         + x(4)*y(1)*z(3) - x(4)*y(3)*z(1) ...
                         - x(2)*y(3)*z(4) + x(2)*y(4)*z(3) ...
                         + x(3)*y(2)*z(4) - x(3)*y(4)*z(2) ...
                         - x(4)*y(2)*z(3) + x(4)*y(3)*z(2) ) /6;
    end
    % Assemble M, K, and s
    for e=1:N_e
        Nodes    = Elements(e,:);
        x        = Points(Nodes,1);
        y        = Points(Nodes,2);
        z        = Points(Nodes,3);
        G        = [((y(4)-y(2))*(z(3)-z(2))-(y(3)-y(2))*(z(4)-z(2))),... 
                    ((y(3)-y(1))*(z(4)-z(3))-(y(3)-y(4))*(z(1)-z(3))),...
                    ((y(2)-y(4))*(z(1)-z(4))-(y(1)-y(4))*(z(2)-z(4))),...
                    ((y(1)-y(3))*(z(2)-z(1))-(y(1)-y(2))*(z(3)-z(1)));...
                    ((x(3)-x(2))*(z(4)-z(2))-(x(4)-x(2))*(z(3)-z(2))),...
                    ((x(4)-x(3))*(z(3)-z(1))-(x(1)-x(3))*(z(3)-z(4))),...
                    ((x(1)-x(4))*(z(2)-z(4))-(x(2)-x(4))*(z(1)-z(4))),...
                    ((x(2)-x(1))*(z(1)-z(3))-(x(3)-x(1))*(z(1)-z(2)));...
                    ((x(4)-x(2))*(y(3)-y(2))-(x(3)-x(2))*(y(4)-y(2))),...
                    ((x(3)-x(1))*(y(4)-y(3))-(x(3)-x(4))*(y(1)-y(3))),...
                    ((x(2)-x(4))*(y(1)-y(4))-(x(1)-x(4))*(y(2)-y(4))),...
                    ((x(1)-x(3))*(y(2)-y(1))-(x(1)-x(2))*(y(3)-y(1)))];
        for p=1:4
            m           = Nodes(p);
            Gp          = [G(1,p), G(2,p), G(3,p)];
            for q=1:4
                n    	= Nodes(q);
            	Gq      = [G(1,q), G(2,q), G(3,q)];
                M(m,n)  = M(m,n) + rho*C*M_e(p,q)*Omega(e)/20;
                K(m,n)	= K(m,n)  - k*dot(Gp,Gq)/(36*Omega(e));
            end 
        end
    end
    
    % Apply boundary conditions
    for b=1:N_b
        if     strcmp(Boundaries(b).type, 'neumann')
            for f=1:Boundaries(b).N;
                Nodes   = Faces(Boundaries(b).indices(f),:);
                for p=1:3
                    m     = Nodes(p);
                    s(m)  = s(m)  + s_e(p)*Qcpu...
                            *Gamma(Boundaries(b).indices(f))/3;
                end
            end
        elseif strcmp(Boundaries(b).type, 'robin')
             for f=1:Boundaries(b).N;
                 Nodes   = Faces(Boundaries(b).indices(f),:);
                 for p=1:3
                     m       = Nodes(p);
                     for q=1:3
                         n       = Nodes(q);
                         K(m,n)  = K(m,n) - h*Gamma(Boundaries(b).indices(f))*K_f(p,q)/12;
                     end
                     s(m)  	= s(m) + h*Tair*s_e(p)*Gamma(Boundaries(b).indices(f))/3;
                  end
             end
        end        
    end
end
 
