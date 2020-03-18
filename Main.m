clc;
clear;
%% Assignment 3 ELEC 4700 Monte-Carlo/Finite Difference Method
% Adam Heffernan 100977570
% Completed March 17th 2020

%%  
% This assignment combines the finite difference and monte carlo methods to
% simulate a silicon wafer under the influence of an electric field and
% with bottle-necks. As shown in both parts 1 and 3 there is no force on
% the particles in the Y direction. This is because there is no potential
% in the Y-direction. Therefore in this direction the particles rely solely
% on scattering to influence there behaviour. However in the X direction we
% see that the particles are influenced by both scattering and the
% magnitude of the electric field across the semi-conductor. The force on
% each electron is extremely small, this is to be expected as the electric
% field we are applying is quite small and therefore should not cause a
% major force on the particles in question. The current reaches a steady
% state at around 10^-12 seconds. 
%
%
%
%% Part 1
Assignment3_Part1;
%%
% Figure 1 shows the trajectories of a couple of particles in the particle
% population that were plotted for the number of time steps specified. The
% particles are being shifted towards the right and that is because there
% is a voltage applied to the left hand side of the semiconductor. This is
% creating an electric field that is forcing the particle to the right.
% Figure 2 shows the semiconductor temperature and the drift current of
% that semi-conductor over the total time for our simulation. Temperature
% and drift current varry until they reach a steady state around 1x10^-12
% seconds. Figure 3 plots the electron density with no bottle neck present
% in the simulation, Figure 4 plots the temperature density map showing the
% temperatures at various locations on the silicon crystal. 
%% 
%% Part 2
Assignment3_Part2;
%%
% Figure 5 shows the conductivity map of the mesh grid with bottle-necks as
% in assignment 2. Figure 6 shows the potential voltage over the contour of
% the silicon crystal.The potential is highest around the boxes as there is
% very high resistance/low conductivity. Therefore according to ohm's law
% this region should have the highest potential. Figure 7 is showing the
% electric field of the mesh grid when a potential of 0.1V is applied to
% the left side of the silicon wafer. All the electric field lines point
% away from the higher potential and therefore all the electrons that
% experience this path will be forced to the right due to the magnitude of
% the electric field. 
%%
%% Part 3 
Assignment3_Part3;
%%
% Figure 8 is showing the trajectories for 10 of 30000 particles in the
% simulation with the bottle neck box included. Figure 9 models the drift
% current vs time and the temperature of the semi-conductor over the same
% time frame. Figure 10 shown above is showing the density of electrons in 
% the mesh grid. These electrons are bunching up due to the electric field 
% acting on them and the fact that the boxes are acting as a bottle neck. 
% Figure 11 shows the bottle neck vs temperature in the silicon crystal.
%






