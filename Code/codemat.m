clear;



%% Physical parameters

% Domain radius (set to 1 by choosing units of length)

R = 1;

% Number of particles

N = 40;

% Swimming speed (set to 1 by choosing units of time)

v0 = 1;

% Particle radius

a = 0.05;

% Squirmer parameter

beta = -7.5;

% Translational diffusivity

D = 0;

% Rotational diffusivity

Do = 0;



%% Numerical parameters

% Time step

dt = 1e-4;

% Cut-off for -log(epsilon)

lnEps_cr = 5;

% Amplitude of steric interactions

Es = 1;

% Simulation time

Tmax = 100;

% Periodicity of outputs

dt_out = 0.01;



%% Some fixed parameters

% Distance of orientational interactions

do = 3*a;

% Amplitude of orientational interactions

Eo = (3/10)*v0/a;

% Distance of steric interactions

ds = 2^(7/6)*a;



fprintf('Effective surface density: %g\n', N*a^2/R^2)



%% Initial data

X = zeros(1,N);

Y = zeros(1,N);

for i = 1:N

    redo = true;

    while redo

        % Generate uniform positions in the disk by rejection

        X(i) = 2*(R-a)*(rand()-0.5);

        Y(i) = 2*(R-a)*(rand()-0.5);

        if X(i)^2 + Y(i)^2<(R-a)^2

            % Check if the distance to any previous point is less than r

            redo = any((X(1:i-1)-X(i)).^2+(Y(1:i-1)-Y(i)).^2<ds^2);

        else

            redo = true;

        end

    end

end

% Checks

dd = 2*R;

for i = 1:N

    for j=1:(i-1)

        dd = min(dd,sqrt((X(i)-X(j))^2+(Y(i)-Y(j))^2));

    end

end

fprintf('Minimal distance: %g \t limit: %g\n', dd, ds)

% Generate random orientations

Theta = 2*pi*rand(1,N);



%% Plots

figure(1), clf

thpl = linspace(0, 2*pi, 100);

plot(R*cos(thpl),R*sin(thpl),'k','LineWidth',2)

hold on

hs = scatter(X, Y, 'filled');

hq = quiver(X,Y,cos(Theta),sin(Theta),'Color','b','AutoScaleFactor',0.4,...

    'LineWidth',2);

axis equal

xlim(1.1*[-R R]);

ylim(1.1*[-R R]);

set(gca,'Position',[0 0 1 1])

currentunits = get(gca,'Units');

set(gca, 'Units', 'Points');

axpos = get(gca,'Position');

set(gca, 'Units', currentunits);

markerWidth = 2*a/diff(xlim)*axpos(3);

set(hs, 'SizeData', markerWidth^2,'MarkerFaceColor','b')

set(gca, 'Visible', 'off')

strng = ['$\begin{array}{l} a =' sprintf('%g', a) ',\ N=' sprintf('%d',N)...

    '\\[5pt] \beta=' sprintf('%g',beta) '\end{array}$'];

text(-1.3,-0.9,strng,'Interpreter','latex','FontSize',24)

%% Time loop

tout = dt_out;

Stat = [];

for t=dt:dt:Tmax

    % Compute all pairwise distances between particles

    DX = X'-X; % matrix whose element (j,i) is X(j)-X(i)

    DY = Y'-Y; % matrix whose element (j,i) is Y(j)-Y(i)

    Rij = sqrt(DX.^2 + DY.^2);

    Rij(1:(N+1):end) = inf;

    

    % Compute steric forces between particles from Weeks-Chandler-Andersen potential

    [I,J] = find(Rij < ds);

    K = sub2ind([N,N],I,J);

    Tmp = -12*(Es/a)*(2*(2*a./Rij(K)).^13-(2*a./Rij(K)).^7)./Rij(K);

    Fs_x = sparse(I,J,Tmp.*DX(K),N,N);

    Fs_y = sparse(I,J,Tmp.*DY(K),N,N);



    % Compute steric forces between particles and wall from

    % Weeks-Chandler-Andersen potential

    Rpart = sqrt(X.^2+Y.^2);

    Fs_pw = zeros(2,N);

    K = find((R-Rpart) < 2^(1/6)*a);

    Tmp = -24*(Es/a)*(2*(a./(R-Rpart(K))).^13-(a./(R-Rpart(K))).^7)./Rpart(K);

    Fs_pw(1,K) = Tmp.*X(K);

    Fs_pw(2,K) = Tmp.*Y(K);



    % Compute the torque Gamma_ij exerted on the i-th particle

    % by its own flow, when modified by a nearby j-th particle

    [I,J] = find(Rij < 3*a);

    Val = zeros(size(I));

    for ip = 1:length(I)

        ex = DX(I(ip),J(ip))./Rij(I(ip),J(ip));

        ey = DY(I(ip),J(ip))./Rij(I(ip),J(ip));

        lnEps = min(lnEps_cr,-log(Rij(I(ip),J(ip))/a-2));

        Val(ip) = Eo*(1+beta*(cos(Theta(I(ip))).*ex+sin(Theta(I(ip))).*ey)) ...

        .*lnEps.*(ex.*sin(Theta(I(ip))) - ey.*cos(Theta(I(ip))));

    end

    Gamma_ij = sparse(I,J,Val,N,N);



    % Compute torque exerted on particles by the wall

    Gamma_w = zeros(1,N);

    K = find((R-Rpart) < 2*a);

    lnEps = min(lnEps_cr,-log((R-Rpart(K))/a-1));

    ex = X(K)/Rpart(K);

    ey = Y(K)/Rpart(K);

    Gamma_w(K) = 2*Eo*(1+beta*(cos(Theta(K)).*ex+sin(Theta(K)).*ey))...

        .*lnEps.*(sin(Theta(K)).*ex-cos(Theta(K)).*ey);



    % Evolution of positions (Euler step)

    X = X + dt*(v0*cos(Theta)+sum(Fs_x,1)+Fs_pw(1,:)) +sqrt(2*dt*D)*randn(1,N);

    Y = Y + dt*(v0*sin(Theta)+sum(Fs_y,1)+Fs_pw(2,:)) +sqrt(2*dt*D)*randn(1,N);

    % Evolution of orientations

    Theta = Theta + dt*(sum(Gamma_ij+0.25*Gamma_ij',1)+Gamma_w) +sqrt(2*dt*Do)*randn(1,N);



    % Reflective boundary conditions at r = R

    R2 = X.^2+Y.^2;

    K = find(R2>(R-a)^2);

    Phi_part = atan2(Y(K),X(K));

    R_part = sqrt(R2(K));

    X(K) = (2*(R-a)-R_part).*cos(Phi_part);

    Y(K) = (2*(R-a)-R_part).*sin(Phi_part);

    

    if t>=tout

        % Plots

        set(hs,'XData',X,'YData',Y);

        set(hq,'XData',X,'YData',Y,'UData',cos(Theta),'VData',sin(Theta));

        xlim(1.1*[-R R]);

        ylim(1.1*[-R R]);

        pause(0.01);

        % Statistics

        Rpart = sqrt(X.^2+Y.^2);

        DT = delaunayTriangulation(X', Y');

        Psi6 = zeros(1,N);

        for ip=1:N

            attachedTri = DT.vertexAttachments(ip);

            attachedPnt = setdiff(unique(DT.ConnectivityList(cell2mat(attachedTri), :)),ip);

            Psi6(ip) = abs(mean((exp(1i*6*atan2(Y(attachedPnt)-Y(ip),X(attachedPnt)-X(ip))))));

        end

        Stat = [Stat; t mean(Rpart) mean((X.*cos(Theta)+Y.*sin(Theta))./Rpart) mean((X.*sin(Theta)-Y.*cos(Theta))./Rpart) mean((X.*sin(Theta)-Y.*cos(Theta))./Rpart.^2) mean(Psi6)];

        tout = tout+dt_out;

    end

end