% henriettes change
clear;
clc
close all
figure('units','normalized','outerposition',[0 0 1 0.6])
colormap jet


system('mkdir pics')
n = 100;          % number of penguin speed
r = 0.5;         % radius of penguin
temp_init = 1;   % temperture of penguin
T0 = 1;

L = 20;         % domain size
noise = 0.01;   % noise amplitude
dt = 0.1;       % time step size

v0 = 0.1;       % penguin speed
a = 1;          % weighting of prior orientation on heading
f = 10000;         % weighting of soft repulsion on orientation
t = 0.1;        % weighting of temperature gradient on orientation
k = 0.5;        % temp decay
tau = 0.01;     % relaxation constant for angle
s = 1;
tau1 = 0.1;     % temperature relaxation

%theta = [0; pi];
N_temp = zeros(n,2);        % orientation vector of penguins
orient = 2*pi*rand(n, 1);   % orientation of penguins
alpha = zeros(n, 1);        % orientation of penguins
N = [cos(orient), sin(orient)];
T = temp_init*ones(n, 1);   % temperature of penguin
sgn = ones(n, 1);
F = zeros(n, 2);            % orientation vector of penguins
H = zeros(n, 2);            % temperature gradient for each penguin to follow
S = zeros(n, 2);            % visual for each penguin to follow
R = zeros(n, n);             % array for penguin penguin distance;
%pos = [6,8; 9,8];
%pos = [L/2, L/2];
pos = 2*L*rand(n,2) - 0.5*L; % position of penguins
LL = 10*L;
Temp = zeros(LL, LL);
bg_temp = 1;

% plot(pos(:,1), pos(:,2),'.')
% hold on
% quiver(pos(:,1), pos(:,2), N(:,1), N(:,2),.2)
% axis([0 L 0 L])

%% Initial equilibration (huddling) phase
for t = 1:500
    
    if(mod(t,100) == 0)
        disp(t)
        % plot
        subplot(1, 2, 1)
        plot(pos(:,1), pos(:,2), 'o', 'Markersize', 11, 'color','k');
        title(['time = ', num2str(t)])
        %axis([L/4 3*L/4 L/4 3*L/4])
        axis([0 L 0 L])
        drawnow
    end
    
    % penguin-penguin distance
    for i = 1:n
        R(:, i) = sqrt((pos(:, 1) - pos(i, 1)).^2 + (pos(:, 2) - pos(i, 2)).^2);
        R(i, i) = 4*r;
    end
    
    % soft core repulsion
    for i = 1:n
        F(i,:) = [0, 0];
        rr = R(:,i);
        I = find(rr < 2*r);
        if(~isempty(I))
            for j = I'
                F(i,:) = F(i,:) - (2*r - rr(j))^(1/2)*(pos(j,:) - pos(i,:))/rr(j);
            end
        end
    end
    
    % huddling
    for i = 1:n
        S(i,:) = [0, 0];
        rr = R(:,i);
        I = 1:n;
        I(I == i) = [];
        for j = I
            pij = pos(j, :) - pos(i, :);
            S(i,:) = S(i,:) + pij/norm(pij);
        end
    end
    
    % recalulcate orientation of penguins
    N_temp = S + f*F;
    
    for i = 1:n
        N_temp(i,:) = N_temp(i,:)/norm(N_temp(i,:));
    end
    N_temp(isnan(N_temp)) = 0;
    
    % update position
    pos = pos + dt*v0*10*N_temp + 0.01*randn(size(pos));
    
end


time = 0;
go = 1;
while(go == 1)
    
    time = time + 1;
    
    % penguin-penguin distance
    for i = 1:n
        R(:, i) = sqrt((pos(:, 1) - pos(i, 1)).^2 + (pos(:, 2) - pos(i, 2)).^2);
        R(i, i ) = 4*r;
    end
    
    % soft repulsion term
    for i = 1:n
        F(i,:) = [0, 0];
        rr = R(:,i);
        I = find(rr < 2*r);
        if(~isempty(I))
            for j = I'
                F(i,:) = F(i,:) - (2*r - rr(j))^(1/2)*(pos(j,:) - pos(i,:))/rr(j);
            end
        end
    end
    
    % temperature gradient term
    for i = 1:n
        H(i,:) = [0, 0];
        rr = R(:,i);
        I = 1:n;
        I(I == i) = [];
        for j = I
            H(i,:) = H(i,:) + sgn(i)*T(j)*exp(-k*rr(j)^2)*(pos(j,:) - pos(i,:))/rr(j);
        end
    end
    
    % update temperature of penguin
    for i = 1:n
        u = H(i,:);
        T(i) = T(i) - dt*(tau1*(T(i) - norm(u)));
        if(T(i) > T0)
            sgn(i) = -1;
        else
            sgn(i) = 1;
        end
    end
    
    % background temperature
    if(bg_temp == 1)
        for i = 1:LL
            for j = 1:LL
                Temp(j, i) = 0;
                x = i/10;
                y = j/10;
                I = 1:n;
                for jj = I
                    Temp(j, i) = Temp(j, i) + T(jj)*exp(-k*((pos(jj,1) - x)^2 + (pos(jj,2) - y)^2));
                end
            end
        end
    end
    
    % recalulcate orientation of penguins
    N_temp = a*N + f*F + t*H;
    
    for i = 1:n
        N_temp(i,:) = N_temp(i,:)/norm(N_temp(i,:));
    end
    N_temp(isnan(N_temp)) = 0;
    
    % update position
    pos = pos + dt*v0*N_temp;
    
    % relax orientation
    for i = 1:n
        alpha(i) = acos(min(1, dot(N_temp(i,:), N(i,:))/(norm(N_temp(i,:))*norm(N(i,:))))); % angle from heading to inst. orientation
        orient(i) = orient(i) + dt*(tau*(-alpha(i)) + noise*rand);                          % add noise to orientation
        N(i,:) = [cos(orient(i)), sin(orient(i))];                                          % position of penguins
    end
    
    % plot
    if(mod(time, 1) == 0)
        subplot(1, 2, 1)
        plot(pos(:,1), pos(:,2), 'o', 'Markersize', 11, 'color','k');
        title(['time = ', num2str(time)])
        axis([0 L 0 L])
        
        if(bg_temp == 1)
            subplot(1, 2, 2)
            imagesc(flipud(Temp));
            caxis([0 T0])
            colorbar
        end
        drawnow
    end
    
    %print(['pics/penguins',num2str(time),'.png'],'-r200','-dpng','-noui');
    %     pause
    
end

