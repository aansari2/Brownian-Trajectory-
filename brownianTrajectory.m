gas.p = 0.1; % Air Pressure Pa (soft vacuum)
gas.T = 300; % Air Temperature Kelvin
gas.mu = 1.96e-5; % gas. viscosity
gas.R = 296.8; % specific gas. constant of air
gas.u = [0;0;0]; % still air
gas.mass = 28*1.66054e-27; % air molecule mass

particle.diameter = 100*1e-6; % 100 nm
particle.density = 1050; % kg/m^3;
particle.volume = 4/3*pi*(particle.diameter/2)^3; % m^3
particle.mass = particle.density*particle.volume; % kg
particle.T = 300;

initialPosition = [0;0;0];
initialVelocity = [0;0;-1];
tspan = [0 1]; %timespan
dt = 0.001;

Ni = diff(tspan)/dt; %Number of iterations
t = tspan(1);

Nt = 1000; % monte carlo simulation count
data = zeros(2,Nt);
for k = 1:Nt 
    q = zeros(6,Ni);
    q(1:3,1) = initialPosition;
    q(4:6,1) = initialVelocity;
    kb = 1.380649 * 10^-23; % m2 kg s-2 K-1 boltzmanns constant

    for i = 1:Ni 
        u = q(4:6,i); % particle. velocity
        h = gas.mass/(2*kb*gas.T);
        hprime = particle.mass/(2*kb*particle.T);
        %[F_drag,F_brwn] = microscopicForce(h,hprime,gas.p,particle.diameter/2,gas.mass,u);
        [F_drag,F_brwn] = stokesForce(u,gas,particle);
        a_g = [0,0,-9.8]';
        dqdt = [u;...
            a_g + (F_drag + randn(3,1)*F_brwn/sqrt(dt))/particle.mass];
        q(:,i+1) = q(:,i) + dt*dqdt;
        t = t + dt;
    end
    data(:,k) = q(1:2,end);
end

averageDeviationDiameter = 2*vecnorm(std(data, [], 2));

plot(data(1,:), data(2,:),'.')
hold on
plot(averageDeviationDiameter*exp(2i*pi*linspace(0,1)))
hold off
axis equal


title(sprintf('trajectory destinations x-y coordinates. Deviation diameter: %2.2e m',...
    averageDeviationDiameter))

xlabel('x (m)'); ylabel('y (m)')


function [F_drag,F_brwn] = stokesForce(Uvec,gas,particle)

kb = 1.380649 * 10^-23; % m2 kg s-2 K-1 boltzmanns constant

lambda = gas.mu./gas.p.*sqrt(pi*gas.R*gas.T/2); % mean free path
Knp = lambda./particle.diameter; % knudsen number of gas.

c1 = 1.2310; c2= 0.4695; c3= -1.1783;
Cc = 1 + Knp.*(c1+c2*exp(c3./Knp)); % cunningham slip factor
F_drag = (3*pi*gas.mu*particle.diameter.*(Uvec-gas.u))./Cc; % Drag Force

S0 = 216*gas.mu*kb*gas.T./... % brownian diffusion coeff
    (pi^2*particle.diameter.^5*particle.density^2.*Cc);

F_brwn = particle.mass.*sqrt(pi*S0); % Stochastic Brownian Force 

end
