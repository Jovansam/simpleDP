tic;

clear all;
close all;

% Declare model
model = Model;

% Initialise environment
model.envInit();

return;

useEuler = 1;
model.eulerInit(useEuler);

% Solve model
model.solve();

% Plot numerical and analytical policy functions
ixt = 1;
model.plotPolicyC(ixt);
model.plotPolicyA1(ixt);
model.plotValue(ixt);
model.plotMargUtil(ixt);

% Simulate model
%%%%%%%%%%%%%%%%%% HERE %%%%%%%%%%%%%%%%%%%%%%%
% Need to add routine to simulate under uncertainty that also works when
% there is no uncertainty. Ideally then get rid of routine that simulates
% under no uncertainty
model.simulate();

% Draw some plots of the simulations
model.plotCons();
model.plotAssets();

% Clear model
%model.clear();


toc;

% Still to do:
    % Why is the solution based on maximising the value function so bad? It
    % must be because the value function approximation is bad

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
    
    % If the assets grid is set so that the lowest each year is the amount
    % of assets that will guarantee minimal consumption in every future
    % year and there is no income, then an individual on that lower bound
    % effectively has no choice over consumption: he must choose minimal
    % consumption. Any more than that would violate minimal assets next
    % period, any less than that would violate the minimal consumption
    % constraint this period
    
    % It is possible to store the policy function in the model solution to
    % remove the need to calculate the optimal policy again in the
    % simulation
    
    % The positioning of the assets grid points can matter a LOT when the
    % solution is not exact
    
    % Suppose you're maximising the value function (i.e. you're not using
    % the Euler equation)
        % In the solution, all you need is the one-period-ahead EV
        % In the simulation, you either need policyA1 or you need EV to 
        % calculate policyA1
        % This means you can either store EV or policyA1 in the solution
        % Is one better than the other? Were we to have a discrete decision
        % how would interpolating the policy function work?
    
    % Useful validation checks:
        % Solution vs exact solution
        % Value function solution vs Euler solution
        % How solution changes as number of grid points increases
        % Value and marginal utility are increasing in the right things
        % Value in last period equals utility (solution and simulation)

% Questions

    % What is checkSDF (check stochastic discount factors)?
    % How would I go about writing unit tests for this model?
