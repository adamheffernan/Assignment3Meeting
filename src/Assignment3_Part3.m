function Assignment3_Part3

len = 200e-9;
width = 100e-9;
electron_charge = -1.60217662e-19; % Charge on electron
voltage_x = 1; % Voltage applied across L (x = 0 is positive)
voltage_y = 0; % Voltage applied across W (y = 0 is positive)
electron_density = 1e15*100^2; % Concentration of electrons in 1/m^2
electron_mass = 9.10938356e-31;
effective_mass = 0.26*electron_mass;
temperature = 300;
k = 1.38064852e-23;
thermal_voltage = sqrt(2*k*temperature/effective_mass);
particle_population = 30000;
particle_plot = 10;
dt = width/thermal_voltage/100;
iter = 200;
scattering_probability = 1 - exp(-dt/0.2e-12);
thermal_distribution = makedist('Normal', 'mu', 0, 'sigma', sqrt(k*temperature/effective_mass));
total_scatters = 0;
top_specular = 0;
bottom_specular = 0;



Ex = voltage_x/len;
Ey = voltage_y/width;
fprintf('Electric Field experienced in the X direction Part 3 %f MV/m\n',Ex/10^6);


Fx = electron_charge*Ex;
Fy = electron_charge*Ey;
fprintf('Force experienced in the X direction Part 3 %f fN\n',Fx/10^-15);


%%
% For one time step, this increases the speed in each direction by

dvx = Fx*dt/effective_mass;
dvy = Fy*dt/effective_mass;
dvx = dvx.*ones(particle_population,1);
dvy = dvy.*ones(particle_population,1);


pos_velo = zeros(particle_population, 4);
trajectories = zeros(iter, particle_plot*2);
temperature = zeros(iter,1);
current_density = zeros(iter,2); % Current density as [Jx Jy] rows


% Generate an initial population
%Creating Box Parameters for boxes shown in plots
x1_box1=80;
x2_box1=120;
y1_box1=0;
y2_box1=40;
x1_box2=80;
x2_box2=120;
y1_box2=60;
y2_box2=100;
%Creating Box Parameters
boxes = 1e-9.*[x1_box1 x2_box1 y1_box1 y2_box1; x1_box2 x2_box2 y1_box2 y2_box2];
%Boxes are specular or diffusive
boxes_specular = [0 1];
for j = 1:particle_population
    %Creating a position, velocity vector in the x and y directions for
    %the number of particles in the simulation.
   
    pos_velo(j,:) = [len*rand width*rand random(thermal_distribution) random(thermal_distribution)];
    while(in_box(pos_velo(j,1:2), boxes))
        %Checking for Particles in box, when particles are in
        %box we will reset there x and y position to a random
        %value
        pos_velo(j,1:2) = [len*rand width*rand];
    end
end
 x_plot1 = [x1_box1, x2_box1, x2_box1, x1_box1, x1_box1];
    y_plot1 = [y1_box1, y1_box1, y2_box1, y2_box1, y1_box1];
    x_plot2 = [x1_box2, x2_box2, x2_box2, x1_box2, x1_box2];
    y_plot2 = [y1_box2, y1_box2, y2_box2, y2_box2, y1_box2];

figure(8);
plot([],[]);
%Creating the box in the plot
%Plotting the Boxes
plot(x_plot1,y_plot1,'-','Color',[0 0 0])
hold on
plot(x_plot2,y_plot2,'-','Color',[0 0 0])
axis([0 len/1e-9 0 width/1e-9]);
title(sprintf('Trajectories for %d of %d Electrons',...
    particle_plot, particle_population));
xlabel('x (nm)');
ylabel('y (nm)');

figure(9);
subplot(2,1,1);
temperature_plot = animatedline;
temperature_plot.Color = 'blue';
title('Temperature');
xlabel('Time');
ylabel('Temperature (K)');
grid on;

figure(9);
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
    for j=1:particle_population
        %Checking for each particle in the population to see if it has hit
        %a boundary of the box we created above. 
        box_num = in_box(pos_velo(j,1:2), boxes);
        x = 0;
        updated_x = 0;
        y = 0;
        updated_y = 0;
        while(box_num ~= 0)
           %If the electron come back that it has hit a boundary then we
           %need to see what boundary it hit in order to predict its next
           %move.
            
            if(pos_velo(j,3) > 0)
                x = pos_velo(j,1) - boxes(box_num,1);
                updated_x = boxes(box_num,1);
            else
                x = boxes(box_num,2) - pos_velo(j,1);
                updated_x = boxes(box_num,2);
            end
            
            if(pos_velo(j,4) > 0)
                y = pos_velo(j,2) - boxes(box_num, 3);
                updated_y = boxes(box_num, 3);
            else
                y = boxes(box_num, 4) - pos_velo(j,2);
                updated_y = boxes(box_num, 4);
            end
            
            if(x < y)
                pos_velo(j,1) = updated_x;
                if(~boxes_specular(box_num))
                    %Diffusive Boundaries in the x direction
                    change = -sign(pos_velo(j,3));
                    v = sqrt(pos_velo(j,3).^2 + pos_velo(j,4).^2);
                    angle = rand()*2*pi;
                    pos_velo(j,3) = change.*abs(v.*cos(angle));
                    pos_velo(j,4) = v.*sin(angle);
                else
                    % Specular Boundaries in the x direction
                    pos_velo(j,3) = -pos_velo(j,3);
                end
            else
                pos_velo(j,2) = updated_y;
                if(~boxes_specular(box_num))
                    %Diffusive Boundaries in the y direction
                    change = -sign(pos_velo(j,4));
                    v = sqrt(pos_velo(j,3).^2 + pos_velo(j,4).^2);
                    angle = rand()*2*pi;
                    pos_velo(j,3) = v.*cos(angle);
                    pos_velo(j,4) = change.*abs(v.*sin(angle));
                else
                    % Specular Boundaries in the y direction
                    pos_velo(j,4) = -pos_velo(j,4);
                end
            end
            %Getting the box that the particle is hitting
            box_num = in_box(pos_velo(j,1:2), boxes);
        end
    end
    % Record the temperature
    temperature(i) = (sum(pos_velo(:,3).^2) + sum(pos_velo(:,4).^2))*effective_mass/k/2/particle_population;
    
    % Record positions and velocities for subset of particles that will be graphed
    for j=1:particle_plot
        trajectories(i, (2*j):(2*j+1)) = pos_velo(j, 1:2);
    end
    

    current_density(i, 1) = electron_charge.*mean(pos_velo(:,3));
    current_density(i, 2) = electron_charge.*mean(pos_velo(:,4));

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
fprintf('Average Electron Velocities in Silicon Crystal for part 3 are %f km/s \n',velocity_average/10^3)
%Output the average time between colisions
fprintf('Average time for part 3 is %f ns\n',t_mn/10^-9);
%Output the mean free path
fprintf('Mean free Path for part 3 is %f mm\n\n',mfp/10^-3);

figure(8);
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
[x y] = meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
f=f./sum(f(:));
figure(10);
electron_density = conv2(electron_density,f,'same');
electron_density = electron_density/(width./size(electron_density,1)*len./size(electron_density,2));
surf(conv2(electron_density,f,'same'));
colorbar
title('Electron Density Map');
xlabel('x (nm)');
ylabel('y (nm)');

%%
% The temperature map is created using a similar procudure. The electron
% velocities are put into bins over space to calculate the temperature at
% different points:
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

%%
% Like with the density map, perform some smoothing:

N = 20;
sigma = 1.5;
[x, y] = meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
f=f./sum(f(:));
figure(11);
surf(conv2(temp,f,'same'));
colorbar
title('Temperature Density Map');
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
