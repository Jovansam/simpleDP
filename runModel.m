tic;

clear all;
close all;

% Declare model
model = Model;

% Initialise environment
model.envInit();

% Solve model
model.solve();

% Find analytical (exact) policy functions
model.solveExact();

% Plot numerical and analytical policy functions
model.plotPolicyC(79);
model.plotPolicyA1(79);
model.plotValue(79);

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
% Check whether the answer is right
    % consumption jumps up at the end of life by more than it does in
    % Cormac's code (why does it jump up in either version?)
% Add error checking

% Discoveries
    % Linear interpolation for the objective function is very bad - it
    % gives extremely poor consumption profiles (all consumption right at
    % end of life)
    % You can write the objective function to be evaluated at any value of
    % the state. It is the job of the objective function to deal with
    % interpolation etc