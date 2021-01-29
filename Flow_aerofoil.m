
%%  D
clc
clear all
% A) Parameters
R = 1;                  % Radius of cylinder
x_c = -0.076;           % Center shift => x
y_c = 0;                % Center shift => y
shift = x_c + i*y_c;    % Center shift in imaginary form
lambda = 0.9;         % R - sqrt(x_c^2 - y_c^2)
psi=0.28;               % Given stream value
delta = 0.003;          % Tolerance
 
% B) Aerofoil(Joukowski)
z_aero_points =[];
 
for theta=-pi:pi/100:pi
    z_aero = (R*cos(theta) + x_c)+i*(R*sin(theta)+y_c);
    z_aero_t = z_aero + ((lambda^2)/z_aero);    % Joukowski transform
    z_aero_points = [z_aero_points;z_aero_t];
end

% Plot aerofoil in real-imaginary plane
plot(real(z_aero_points),imag(z_aero_points),'Color','k','LineWidth',2)
grid on
hold on 

% C) Stream
stream = [];
 
% Transforming each and every point of stream function in xy plane to real-imaginary plane
% i) For psi=0 (on the cylinder)
for r=1:0.0005:2
    z_1 = r + shift;
    z_2 = -r + shift;
    
    z_i = [z_1,z_2,z_1,z_2];
    z_i_t= z_i + ((lambda^2)./z_i); % Joukowski transform
    stream = [stream;z_i_t];
end

% ii) For psi~=0(not on cylinder)
for a=0.05:0.25:2.5
    for r=1:0.0005:2
        psi_a = a*psi;
        theta_s = real(asin((psi_a.*r)./(r.^2 - 1))); % Basic angle derived from stream function
        z = r*cos(theta_s) + i*(r*sin(theta_s));
        
        % All possible values of z(4 quadrants)
        z_1 = z + shift;
        z_2 = conj(-z) + shift;
        z_3 = conj(z) + shift;
        z_4 = -z + shift;
                
        % Tolerance(Remove points within cylinder or aerofoil)
        if abs(abs(z_1)-R) <= delta
            z_1 = NaN;
        end
        if abs(abs(z_2)-R) <= delta
            z_2 = NaN;
        end
        if abs(abs(z_3)-R) <= delta
            z_3 = NaN;
        end
        if abs(abs(z_4)-R) <= delta
            z_4 = NaN;
        end
        
        % Joukowski transform of points of stream function
        z_i = [z_1,z_2,z_3,z_4];
        z_i_t= z_i + ((lambda^2)./z_i); 
        
        % Note: Points on each stream function is stored in different rows
        % Thus, we'll have n rows and 4 columns for our variable 'stream'
        stream = [stream;z_i_t];
    end
end
 
% Plot transformed stream function in real-imaginary plane
plot(real(stream),imag(stream),'b')
axis equal
% Labels
xlabel("Real")
ylabel("Imaginary")
title("NACA0012 aerofoil in uniform flow(\Psi = \pm 0.28) (angle of attack,\alpha=0)")
hold off
 


