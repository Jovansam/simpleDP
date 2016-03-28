tic;

clear all;
close all;

% Declare model
model = Model;

% Initialise environment
model.envInit();
useEuler = 1;
model.eulerInit(useEuler);

% Solve model
model.solve();

% Find analytical (exact) policy functions
model.solveExact();

% Plot numerical and analytical policy functions
model.plotPolicyC(1);
model.plotPolicyA1(1);
model.plotValue(1);
model.plotMargUtil(1);

% Simulate model
model.simulate();
model.simulateExact();

% Draw some plots of the simulations
model.plotCons();
model.plotAssets();

% Clear model
%model.clear();


toc;

% Still to do:
    % Add error checking
    % Sort out graphing code (e.g. graph numbering, create graphs for
    % simulated value and marginal utility)

% To ask Cormac

    % 1. Why does consumption jump up at the end of life in the value function
    % solution?
    % 2. Why is the value function more accurately calculated in the value
    % function solution than the Euler solution?
    % 3. What is the reason for the ubA1 - lbA1 < minCons condition in the
    % solveEulerEquation function?

% Discoveries
    
    % Linear interpolation for the value function is very bad - it
    % gives extremely poor consumption profiles (all consumption right at
    % end of life)
    
    % If you solve the Euler equation to find optimal consumption, then the
    % value function is not used. This means that the quality of the value
    % function approximation may not matter
    
    % You can write the objective function to be evaluated at any value of
    % the state. It is the job of the objective function to deal with
    % interpolation etc
    
    % If the assets grid is set so that the lowest each year is the
    % amount of assets that will guarantee minimal consumption in every
    % future year, then an individual on that lower bound effectively has
    % no choice over consumption: he must choose minimal consumption. Any
    % more than that would violate minimal assets next period, any less
    % than that would violate the minimal consumption constraint this
    % period