clear all; close all; clc
% This script aims to compute the Lorentz Force acting on a charged particle
% in a Penning Trap, solving the accompanying kinematics equations to
% ultimately simulate an animation of the particle's motion within the trap.
% All values shown are in SI units.

% =========================================================================
% INPUT PARAMETERS (adjust as necessary)
% =========================================================================

% Proton,  q = 1.602e-19 C    m = 1.672e-27 kg
q = 1.602e-19;   m = 1.672e-27;

% Potential difference between electrodes
V0 = 35.75;

% Defining the dimensions of the Penning Trap (hyperbolic endcaps and ring electrodes, Fig 5 in accompanying write-up)
z0 = 0.01215;   p0 = 0.015;

% Quadrupole potential set up by the electrodes
% V=(V0/(z0^2)+1/2*p0^2).*[-1/2*x^2 -1/2*y^2 z^2];

% Strength of B-field in +z direction (>SQRT((2m/q)V0/z0^2) for stability),
% see accompanying write-up
B = [0 0 0.2];

% Initial velocity
u = [300 400 50];

% Initial displacement, ensure values are within the confines of the trap, z0 and p0
s = [0.001 0.001 0.001];

% Total number of time steps
N = 10000; % (default)

% Time step size, see accompanying write-up
h = abs(0.01*m/(q*norm(B))); % (default) Allows the Reduced Cyclotron Motion 
                             %           to be more pronounced, at the cost
                             %           of long plotting time before the
                             %           superposition of cyclotron orbit, axial
                             %           oscillation and magnetron drift of
                             %           the particle's motion can be
                             %           observed. Increase value 
                             %           by orders of 10 to make the
                             %           Magnetron drift more
                             %           prominent (eg. try 5e-9 and 5e-8).
% =========================================================================
% SOLVING WITH FINITE DIFFERENCE METHOD
% =========================================================================
% Reference: http://www.physics.usyd.edu.au/teach_res/mp/doc/em_vBE.pdf

% Initialising arrays
x  = zeros(N,1);   y  = zeros(N,1);   z  = zeros(N,1);

% Time step 1: n = 1   
x(1) = s(1);   y(1) = s(2);   z(1) = s(3);

% Time step 2: n = 2
x(2) = x(1) + u(1)*h;   y(2) = y(1) + u(2)*h;   z(2) = z(1) + u(3)*h;

%Plotting initial displacement
i=1:2;
plot3(x(i), y(i), z(i), '.r');
xlabel('x'); ylabel('y'); zlabel('z');
drawnow;
grid on;
hold on;

% =========================================================================
% TIME LOOPS for time step 3 onwards
% =========================================================================

% Defining the constant
k = q*norm(B)*h/(2*m);

for n = 2 : N-1
    
    % Displacement
    if abs(x(n)) < p0 && abs(y(n)) < p0 && abs(z(n)) < z0 %checking that the particle is within the confines of the trap
        
        % E-field from the quadrupole potential, see accompanying write-up
        E = -(V0/(z0^2)+1/2*p0^2).*[-x(n) -y(n) 2*z(n)];
        x(n+1) = (2*x(n)+(k^2-1)*x(n-1)+2*k*y(n)-2*k*y(n-1)+(q*E(1)*h^2/m)+k*(q*E(2)*h^2/m))/(1+k^2);
        y(n+1) = (2*y(n)+(k^2-1)*y(n-1)-2*k*x(n)+2*k*x(n-1)+(q*E(2)*h^2/m)-k*(q*E(1)*h^2/m))/(1+k^2);
        z(n+1) = 2*z(n)-z(n-1)+(q*E(3)*h^2/m);
        
        % Plotting the motion
        plot3(x(n+1), y(n+1), z(n+1), '.r');
        xlabel('x'); ylabel('y'); zlabel('z');
        drawnow;
        grid on;
        hold on;
        
        % Saving the animation as a gif
        % Reference: https://ntulearn.ntu.edu.sg/bbcswebdav/pid-1357298-dt-content-rid-5090941_1/courses/17S2-PH2102-LEC/Animations%20in%20MATLAB.pdf?target=blank
        current_frame = getframe(gcf);
        if n == 2
            [mov(:,:,1,n), map] = rgb2ind(current_frame.cdata, 256, 'nodither');
        else
            mov(:,:,1,n) = rgb2ind(current_frame.cdata, map, 'nodither');
        end
        
    else
        disp('Collision of particle with Trap walls!')
        break
    end
    
    imwrite(mov, map, 'PH2102_PenningTrap_animation_jeremylian.gif', 'DelayTime', 0.03, 'LoopCount', inf);
    
end