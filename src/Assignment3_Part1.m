function Assignment3_Part1()
total_scatters = 0;
len = 200e-9;
width = 100e-9;
electron_charge = -1.60217662e-19; % Charge on electron
voltage_x = 0.1; % Voltage applied across L (x = 0 is positive)
voltage_y = 0; % Voltage applied across W (y = 0 is positive)
electron_density = 1e15*100^2; % Concentration of electrons in 1/m^2
electron_mass = 9.10938356e-31;
effective_mass = 0.26*electron_mass;
temperature = 300;
k = 1.38064852e-23;
thermal_velocity = sqrt(2*k*temperature/effective_mass);
mfp = thermal_velocity*0.2e-12; % Mean free path
particle_population = 30000;
particle_plot = 10;
dt = width/thermal_velocity/100;
iter = 200;
scattering_probability = 1 - exp(-dt/0.2e-12);
thermal_distribution = makedist('Normal', 'mu', 0, 'sigma', sqrt(k*temperature/effective_mass));


top_specular = 0;
bottom_specular = 0;


Ex = voltage_x/len;
Ey = voltage_y/width;
fprintf('Electric Field experienced in the X direction Part 1 %f MV/m\n',Ex/10^6);


Fx = electron_charge*Ex;
Fy = electron_charge*Ey;
fprintf('Force experienced in the X direction Part 1 %f fN\n',Fx/10^-15);



dvx = Fx*dt/effective_mass;
dvy = Fy*dt/effective_mass;
dvx = dvx.*ones(particle_population,1);
dvy = dvy.*ones(particle_population,1);

pos_velo = zeros(particle_population, 4);
trajectories = zeros(iter, particle_plot*2);
temperature = zeros(iter,1);
current_density = zeros(iter,2); % Current density as [Jx Jy] rows

for i = 1:particle_population
    pos_velo(i,:) = [len*rand width*rand random(thermal_distribution) random(thermal_distribution)];
end


figure(2);
subplot(2,1,1);
temperature_plot = animatedline;
temperature_plot.Color = 'blue';
title('Temperature');
xlabel('Time');
ylabel('Temperature (K)');
grid on;

figure(2);
subplot(2,1,2);
current_plot = animatedline;
current_plot.Color = 'red';
title('Current Density');
xlabel('Time');
ylabel('Current density (A/m)');
grid on;

%%
% Run through the simulation:
for i = 1:iter
    % Update the velocities
    
    pos_velo(:,3) = pos_velo(:,3) + dvx;
    pos_velo(:,4) = pos_velo(:,4) + dvy;
    
    %Update the positions
    pos_velo(:,1:2) = pos_velo(:,1:2) + dt.*pos_velo(:,3:4);

    j = pos_velo(:,1) > len;
    pos_velo(j,1) = pos_velo(j,1) - len;
    
    j = pos_velo(:,1) < 0;
    pos_velo(j,1) = pos_velo(j,1) + len;
    
    j = pos_velo(:,2) > width;

    if(top_specular)
        pos_velo(j,2) = 2*width - pos_velo(j,2);
        pos_velo(j,4) = -pos_velo(j,4);
    else % Diffusive
        % The electron bounces off at a random angle
        pos_velo(j,2) = width;
        v = sqrt(pos_velo(j,3).^2 + pos_velo(j,4).^2);
        angle = rand([sum(j),1])*2*pi;
        pos_velo(j,3) = v.*cos(angle);
        pos_velo(j,4) = -abs(v.*sin(angle));
    end
    
    j = pos_velo(:,2) < 0;
    
    if(bottom_specular)
        pos_velo(j,2) = -pos_velo(j,2);
        pos_velo(j,4) = -pos_velo(j,4);
    else % Diffusive
        % The electron bounces off at a random angle
        pos_velo(j,2) = 0;
        v = sqrt(pos_velo(j,3).^2 + pos_velo(j,4).^2);
        angle = rand([sum(j),1])*2*pi;
        pos_velo(j,3) = v.*cos(angle);
        pos_velo(j,4) = abs(v.*sin(angle));
    end
    
    % Scatter particles
    j = rand(particle_population, 1) < scattering_probability;
    pos_velo(j,3:4) = random(thermal_distribution, [sum(j),2]);
    total_scatters = total_scatters + count_scatters(i) ;
    
    % Record the temperature
    temperature(i) = (sum(pos_velo(:,3).^2) + sum(pos_velo(:,4).^2))*effective_mass/k/2/particle_population;
    
    % Record positions and velocities for subset of particles that will be graphed
    for j=1:particle_plot
        trajectories(i, (2*j):(2*j+1)) = pos_velo(j, 1:2);
    end
    
    % Calculate and record the current density using the formula from above
    current_density(i, 1) = electron_charge.*mean(pos_velo(:,3));
    current_density(i, 2) = electron_charge.*mean(pos_velo(:,4));
    
    % Plot the temperature and current
    addpoints(temperature_plot, dt.*i, temperature(i));
    addpoints(current_plot, dt.*i, current_density(i,1));

  
end
%Calculate the average time
t_mn = ((dt*iter*particle_population)/total_scatters);
%Calculate the average velocity
velocity_average= sqrt(sum((pos_velo(:,3).^2))/particle_population + sum(pos_velo(:,4).^2)/particle_population);
%Calculate the mean free path
mfp = (t_mn * velocity_average);
%Output the average velocity
fprintf('Average Electron Velocities in Silicon Crystal for part 1 are %f km/s \n',velocity_average/10^3)
%Output the average time between colisions
fprintf('Average time for part 1 is %f ns\n',t_mn/10^-9);
%Output the mean free path
fprintf('Mean free Path for part 1 is %f mm\n\n',mfp/10^-3);
% Show trajectories after the movie is over
figure(1);
title(sprintf('Electron Trajectories for %d of %d Electrons',...
    particle_plot, particle_population));
xlabel('x (nm)');
ylabel('y (nm)');
axis([0 len/1e-9 0 width/1e-9]);
grid on;
hold on;
for i=1:particle_plot
    plot(trajectories(:,i*2)./1e-9, trajectories(:,i*2+1)./1e-9, '.');
end


electron_density = hist3(pos_velo(:,1:2),[200 100])';
N = 20;
sigma = 1.5;
[x, y] = meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
f=f./sum(f(:));
figure(3);
electron_density = conv2(electron_density,f,'same');
surf(conv2(electron_density,f,'same'));
colorbar
title('Electron Density Map');
xlabel('x (nm)');
ylabel('y (nm)');

temp_sum_x = zeros(ceil(len/1e-9),ceil(width/1e-9));
temp_sum_y = zeros(ceil(len/1e-9),ceil(width/1e-9));
temp_num = zeros(ceil(len/1e-9),ceil(width/1e-9));

% Look at velocities of all the particles
for i=1:particle_population
    % Find which "bin" it belongs in:
    x = floor(pos_velo(i,1)/1e-9);
    y = floor(pos_velo(i,2)/1e-9);
    if(x==0)
        x = 1;
    end
    if(y==0)
        y= 1;
    end
    
    % Add its velocity components to the cumulative count:
    temp_sum_y(x,y) = temp_sum_y(x,y) + pos_velo(i,3)^2;
    temp_sum_x(x,y) = temp_sum_x(x,y) + pos_velo(i,4)^2;
    temp_num(x,y) = temp_num(x,y) + 1;
end

%%
% Now, with the velocities added up, calculate the temperatures:
temp = (temp_sum_x + temp_sum_y).*effective_mass./k./2./temp_num;
temp(isnan(temp)) = 0;
temp = temp';

N = 20;
sigma = 1.5;
[x, y] = meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
f=f./sum(f(:));
figure(4);
surf(conv2(temp,f,'same'));
title('Temperature Density Map');
colorbar
xlabel('x (nm)');
ylabel('y (nm)');
end
function count=count_scatters(matrix)
%This function creates a counter which counts the number of colisions of
%particles.
count = 0;
for i=1:size(matrix)
    if matrix(i) == 1
        count=count+1;
    end
end
end