
% Van Der Pol parameters
mu = 1.0;               % see wiki page

% Phase portrait parameters
Nx = 5;                 % Quiver grid extent
dx = 0.5;               % Quiver grid resolution

% Trajectory parameters
tmax = 20;              % Trajectory time-span





% Generate phase portrait
[x, y] = meshgrid(-Nx:dx:Nx, -Nx:dx:Nx);    % Setup grid
xdot = y;                                   % xdot (http://www.math.psu.edu/melvin/phase/newphase.html)
ydot = (1-x.^2).*y-x;                       % ydot (                        "                         )
figure; quiver(x,y,xdot, ydot,'k');         % Draw quiver plot (arrows)

% Simulate trajectories
for n = (1:100)
    y0 = 8*(rand(2,1)-0.5);                 % Initial conditions
    tspan = [0, tmax];                      % Time span
    ode = @(t,y) vanderpoldemo(t,y,mu);     % Matlab comes with a Van Der Pol demo function
    [t,y] = ode45(ode, tspan, y0);          % Ordinary differential equation solver (numerical)
    hold('on'); plot(y(:,1),y(:,2),'k');    % Overlay trajectory onto quiver
    plot(y(1,1),y(1,2),'k*');               % Highlight start point
end;
