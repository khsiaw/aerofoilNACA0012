%% C
clc
clear all
% A) Parameters
R = 1;                  % Radius of cylinder
x_c = -0.076;           % Center shift => x
y_c = 0;                % Center shift => y
shift = x_c + i*y_c;    % Center shift in imaginary form
lambda = 0.9;           % R - sqrt(x_c^2 - y_c^2)
psi=0.28;               % Given stream value
delta = 0.003;          % Tolerance

% Cylinder plot
cylinder = [];
for theta=-pi:pi/100:pi
    cylinder=[cylinder;R*cos(theta),R*sin(theta)];
end
plot(cylinder(:,1),cylinder(:,2),'LineWidth',2,'Color','k')
hold on
grid on
axis equal

% Stream function
xy_coor = [];

% At psi=0
for r=1:0.0005:2
    z = r+i*0;
    xy_coor=[xy_coor;r,r,-r,-r];
end

% At psi~=0
for a=0.05:0.25:2.5
    for r=1:0.0005:2
        psi_a = a*psi;
        
        % Basic angle derived from stream function
        theta_s = real(asin((psi_a.*r)./(r.^2 - 1))); 
        
        z_1 = r*cos(theta_s)+i*(r*sin(theta_s));
        z_2 = conj(z_1);
        z_3 = -z_1;
        z_4 = conj(-z_1);
        
        % Tolerance(Remove points within cylinder or aerofoil)
        if abs(r-R) <= delta
            z_1 = NaN;
            z_2 = NaN;
            z_3 = NaN;
            z_4 = NaN;
        end
        xy_coor = [xy_coor;z_1,z_2,z_3,z_4];       
    end
end

plot(real(xy_coor),imag(xy_coor),'b')
hold off
xlabel("x(m)")
ylabel("y(m)")
title("Uniform flow across cylinder(\Psi = \pm 0.28) with angle of attack,\alpha = 0")


