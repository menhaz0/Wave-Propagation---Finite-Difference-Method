%% Discrete Cutoff
% This first section is a wave propagating through two discrete mediums


% camera stuff
az = 0; % azimuthal angle
el = 30; % elevation

% disctretize oud domain (unit square)
J = 100;
x = linspace(0,1,J+1);
y = linspace(0,1,J+1);

% set up the needed variables for finite difference method
dx = 1/J;
dy = dx;
s = 1/3; % theory says s < 1/sqrt(2) for stability, this feels safe
dt = s*dx;
n1 = 1; % refractive index of first area
n2 = 2; % refractive index of second area
cutoff = floor(J/2); % determines at which point the areas change

% initial data is all zero, I will use a source term
phi = zeros(J+1,J+1);
U = phi;
U_old = phi;




for n = 1:300

    % make matrix for next iteration
    U_new = zeros(J+1,J+1);
    
    for i = 2:J % boundaries set to zero, but not quite
        for j = 2:cutoff % our equation changes between mediums
            U_new(i,j) = (s/n1)^2* ... % this is just s*v, here v = 1/n1
                         ...
                     ... % centered difference laplacian
                         ( U(i+1,j) + U(i,j+1) + U(i-1,j) + U(i,j-1)...
                         - 4*U(i,j)) ...
                         + 2*U(i,j) - U_old(i,j) ...
                         ...
                     ... % contribution from source term
                         + S(j*dx,i*dx,n*dt); 

            % dampen funtion near edges
            U_new(i,j) = U_new(i,j)*D(j*dx,i*dx,J); 
        end
        for j = cutoff:J % now we're in the new medium
            U_new(i,j) = (s/n2)^2* ... % new refractive index, v = 1/n2
                         ( U(i+1,j) + U(i,j+1) + U(i-1,j) + U(i,j-1)...
                         - 4*U(i,j)) ...
                         + 2*U(i,j) - U_old(i,j) ...
                         + S(j*dx,i*dx,n*dt);

            U_new(i,j) = U_new(i,j)*D(j*dx,i*dx,J);
        end
    end


    % plotting time
    hold off

    % I switched my axes, so I had to transpose the matrix
    surf(x, y, U_new', 'EdgeColor','none') 

    axis([0 1 0 1 -0.5 0.5])
    view(az,el); % sets elevation and angle of camera
    az = az - 0.5; % rotates the angle to make it dynamic
    el = el + 0.2;

    % update for next iteration
    U_old = U;
    U = U_new;
    pause(0.05);
end

% this pause is so it doesn't go off to the next section
pause(5);

%% Continuous Change
% what if the refractive index gradually changed? 
% like entering an atmosphere

az = 0;
el = 30;

J = 100;
x = linspace(0,1,J+1);
y = linspace(0,1,J+1);

dx = 1/J;
dy = dx;
s = 1/3;
dt = s*dx;
n1 = 1;
n2 = 10;
cutoff = floor(J/2);

phi = zeros(J+1,J+1);
U = phi;
U_old = phi;




for n = 1:300

    U_new = zeros(J+1,J+1);
    
    for i = 2:J
        for j = 2:J
            % now our refractive index changes for each j, from n1 to n2
            % if I was smart I'd sum over j first, but that looks ugly
            v = 1/(j*dx*(n2-n1) + n1); 

            U_new(i,j) = (s*v)^2* ...
                         ( U(i+1,j) + U(i,j+1) + U(i-1,j) + U(i,j-1)...
                         - 4*U(i,j)) ...
                         + 2*U(i,j) - U_old(i,j) ...
                         + S(j*dx,i*dx,n*dt);
            U_new(i,j) = U_new(i,j)*D(j*dx,i*dx,J);
        end

    end

    hold off
    surf(x, y, U_new', 'EdgeColor','none')
    axis([0 1 0 1 -0.5 0.5])
    view(az,el);
    az = az - 0.5;
    el = el + 0.2;

    U_old = U;
    U = U_new;
    pause(0.05);
end

pause(5);


%% Discrete Cutoff with reflective barrier
% here I place a reflective wall at x = 0.5, 0 < y < 0.5

az = 0;
el = 30;

J = 100;
x = linspace(0,1,J+1);
y = linspace(0,1,J+1);

dx = 1/J;
dy = dx;
s = 1/3;
dt = s*dx;
n1 = 1.5;
n2 = 1;
cutoff = floor(J/2);

phi = zeros(J+1,J+1);
U = phi;
U_old = phi;




for n = 1:300

    U_new = zeros(J+1,J+1);
    
    for i = 2:J
        for j = 2:cutoff
            U_new(i,j) = (s/n1)^2* ...
                         ( U(i+1,j) + U(i,j+1) + U(i-1,j) + U(i,j-1)...
                         - 4*U(i,j)) ...
                         + 2*U(i,j) - U_old(i,j) ...
                         + S(j*dx,i*dx,n*dt);
            U_new(i,j) = U_new(i,j)*D(j*dx,i*dx,J);
        end
        for j = cutoff:J
            U_new(i,j) = (s/n2)^2* ...
                         ( U(i+1,j) + U(i,j+1) + U(i-1,j) + U(i,j-1)...
                         - 4*U(i,j)) ...
                         + 2*U(i,j) - U_old(i,j) ...
                         + S(j*dx,i*dx,n*dt);

            U_new(i,j) = U_new(i,j)*D(j*dx,i*dx,J);
        end
    end

    % a wall cannot have amplitude on it, and making it zero reflects
    for i = 2:floor(J/2)
        U_new(i,cutoff) = 0;
    end

    hold off
    surf(x, y, U_new', 'EdgeColor','none')
    axis([0 1 0 1 -0.5 0.5])
    view(az,el);
    az = az - 0.5;
    el = el + 0.2;

    U_old = U;
    U = U_new;
    pause(0.05);
end

pause(5);



%% Discrete Cutoff with absorbiing barrier
% now I want the wall to absorb the wave

az = 0;
el = 30;

J = 100;
x = linspace(0,1,J+1);
y = linspace(0,1,J+1);

dx = 1/J;
dy = dx;
s = 1/3;
dt = s*dx;
n1 = 1.5;
n2 = 1;
cutoff = floor(J/2);

phi = zeros(J+1,J+1);
U = phi;
U_old = phi;




for n = 1:300

    U_new = zeros(J+1,J+1);
    
    for i = 2:J
        for j = 2:cutoff
            U_new(i,j) = (s/n1)^2* ...
                         ( U(i+1,j) + U(i,j+1) + U(i-1,j) + U(i,j-1)...
                         - 4*U(i,j)) ...
                         + 2*U(i,j) - U_old(i,j) ...
                         + S(j*dx,i*dx,n*dt);
            U_new(i,j) = U_new(i,j)*D(j*dx,i*dx,J);

            % very ugly solution, I basically haphazardly add
            % my dampening B.C. in the middle of the code
            if (i < floor(J/2) )
                U_new(i,j) = U_new(i,j)*(1 - exp(-abs(j*dx-0.5)*J*4));
            end
        end
        for j = cutoff:J
            U_new(i,j) = (s/n2)^2* ...
                         ( U(i+1,j) + U(i,j+1) + U(i-1,j) + U(i,j-1)...
                         - 4*U(i,j)) ...
                         + 2*U(i,j) - U_old(i,j) ...
                         + S(j*dx,i*dx,n*dt);

            U_new(i,j) = U_new(i,j)*D(j*dx,i*dx,J);
        end
    end

    % this isn't needed but ensures there is no amplitude on the wall
    for i = 2:floor(J/2)
        U_new(i,cutoff) = 0;
    end

    hold off
    surf(x, y, U_new', 'EdgeColor','none')
    axis([0 1 0 1 -0.5 0.5])
    view(az,el);
    az = az - 0.5;
    el = el + 0.2;

    U_old = U;
    U = U_new;
    pause(0.05);
end

pause(5);




%% Abusing the absorbing barrier 
% this barrier can be used to make a wave traveling in one direction
% with this I can show total internal reflection

az = 0;
el = 30;

J = 100;
x = linspace(0,1,J+1);
y = linspace(0,1,J+1);

dx = 1/J;
dy = dx;
s = 1/3;
dt = s*dx;
n1 = 1.5;
n2 = 0.8;
cutoff = floor(J/2);

phi = zeros(J+1,J+1);
U = phi;
U_old = phi;




for n = 1:350

    U_new = zeros(J+1,J+1);
    
    for i = 2:J
        for j = 2:cutoff
            U_new(i,j) = (s/n1)^2* ...
                         ( U(i+1,j) + U(i,j+1) + U(i-1,j) + U(i,j-1)...
                         - 4*U(i,j)) ...
                         + 2*U(i,j) - U_old(i,j) ...
                         + S(j*dx,i*dx,n*dt);
            U_new(i,j) = U_new(i,j)*D(j*dx,i*dx,J);

            % shortening the "barrier" and making its dampening reach 
            % farther than it should kills the +y part of the wave for
            % a moment, making the wave hit the cutoff with a wider angle
            if (i < floor(J/3) )
                U_new(i,j) = U_new(i,j)*(1 - exp(-abs(j*dx-0.5)*J/4));
            end
        end
        for j = cutoff:J
            U_new(i,j) = (s/n2)^2* ...
                         ( U(i+1,j) + U(i,j+1) + U(i-1,j) + U(i,j-1)...
                         - 4*U(i,j)) ...
                         + 2*U(i,j) - U_old(i,j) ...
                         + S(j*dx,i*dx,n*dt);

            U_new(i,j) = U_new(i,j)*D(j*dx,i*dx,J);
        end
    end


    for i = 2:floor(J/2)
        U_new(i,cutoff) = 0;
    end

    hold off
    surf(x, y, U_new', 'EdgeColor','none')
    axis([0 1 0 1 -0.5 0.5])
    view(az,el);
    az = az - 0.5;
    el = el + 0.25;

    U_old = U;
    U = U_new;
    pause(0.05);
end

pause(5);

%% Functions
% here are the two functions I always use


% this is my source term, a basic oscillating gaussian
function z = S(x,y,t); 
    z = 0.005*cos(40*t)*exp(-200*((x-0.25).^2 + (y-0.25).^2));
end


% This is my dampening term, 1 everywhere except when near the edges
% where x = 0 and y = 0. This does mean the far walls reflect the wave
function z = D(x,y,J);
    z = 1 - exp(-J*x/10)*exp(-J*y/10);
end