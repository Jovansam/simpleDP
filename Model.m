classdef Model < handle
    properties (Access = public)
        env
        sln
        slnExact
        sims
        simsExact
    end
    methods (Access = private)

        % Utility function
        function utils = utility(obj, cons)
            if (cons <= 0)
                error('Error in utility: consumption <- 0');
            end
            if (obj.env.gamma == 1)
                utils = log(cons);
            else
                utils = (cons.^(1 - obj.env.gamma)) ./ (1 - obj.env.gamma);
            end
        end
        
        % Marginal utility function
        function margUtils = marginalUtility(obj, cons)
            if (obj.env.gamma == 1)
                margUtils = 1 ./ cons;
            else
                margUtils = cons.^(-obj.env.gamma);
            end
        end
        
        % Inverse marginal utility function
        function invMargUtils = inverseMarginalUtility(obj, margUtils)
            if (obj.env.gamma == 1)
                invMargUtils = 1 ./ margUtils;
            else
                invMargUtils = margUtils.^(-1/obj.env.gamma);
            end
        end
        
        % Objective function
        % Rather than passing ixt (time index) as a separate argument, we
        % could pass a state structure combining ixt and A0
        function value = objectiveFunc(obj, ixt, A0, Y, A1)
            cons = A0 + Y - (A1 ./ (1 + obj.env.r));
            VA1 = interp1(obj.env.Agrid(ixt+1, :), obj.sln.value(ixt+1, :), A1, obj.env.interpMethod, 'extrap');
            value = utility(obj, cons) + (obj.env.beta * VA1);
            value = -value;
        end

        % Function for solving the Euler equation
        % Interpolates inverse marginal utility because this is linear
        function eulerDiff = eulerDifference(obj, ixt, A0, Y, A1)
            invMargUtilAtA1 = interp1(obj.env.Agrid(ixt+1, :), inverseMarginalUtility(obj, obj.sln.margUtil(ixt+1, :)), A1, obj.env.interpMethod, 'extrap');
            margUtilAtA1 = marginalUtility(obj, invMargUtilAtA1);
            cons = A0 + Y - (A1 ./ (1 + obj.env.r));
            eulerDiff = marginalUtility(obj, cons) - (obj.env.beta * (1 + obj.env.r) * margUtilAtA1);
        end
        
    end
    methods (Access = public)
        
        % Function to initialise the environment
        function envInit(obj)
            obj.env.tol = 1e-10;
            obj.env.minCons = 1e-5;
            obj.env.T = 80;
            obj.env.tretire = obj.env.T + 1;
            obj.env.r = 0.01;
            obj.env.beta = 1/(1+obj.env.r);
            obj.env.gamma = 1.5;
            obj.env.mu = 0;                     % Log income
            obj.env.startA = 0;
            obj.env.numPointsA = 20;
            obj.env.interpMethod = 'pchip';     % Previously I had 'linear' but this gave very misleading simulations
            obj.env.isUncertainty = 1;
            obj.env.numPointsY = 2;
            obj.env.hcIncome = [0.1, 0.9]';
            obj.env.hcIncPDF = [0.5, 0.5]';
            obj.YgridInit();
            obj.env.borrowingAllowed = 1;
            obj.AgridInit();
            obj.env.solveUsingEuler = 0;
        end
        
        % Function to initialise either to solve using Euler equation or
        % using Value function
        function eulerInit(obj, solveUsingEuler)
            if (solveUsingEuler == 1)
                obj.env.solveUsingEuler = 1;
                obj.env.numpointsA = 10;
                % Note: setting interpMethod to 'linear' means that linear
                % interpolation is used both for linearised marginal
                % utility in the solution of the Euler equation and in
                % finding the value that corresponds to that solution. The
                % latter is a very poor approximation, but it doesn't
                % matter because we don't use the value function in any of
                % the solution
                obj.env.interpMethod = 'linear';
            else
                obj.env.solveUsingEuler = 0;
                obj.env.numpointsA = 20;
                obj.env.interpMethod = 'pchip';
            end
        end
        
        % Function to populate income grid
        function YgridInit(obj)
            if (obj.env.isUncertainty == 0)
                obj.env.Ygrid = repmat(exp(obj.env.mu), obj.env.T, 1);
                obj.env.Ymin = obj.env.Ygrid;
                obj.env.Ymax = obj.env.Ygrid;
                obj.env.Q = 1;                  % Markov transition matrix
            elseif (obj.env.isUncertainty == 1)
                obj.env.Ygrid = repmat(obj.env.hcIncome', obj.env.T, 1);
                obj.env.Ymin = min(obj.env.Ygrid, [], 2);
                obj.env.Ymin = max(obj.env.Ygrid, [], 2);
            end
            if (obj.env.tretire <= 0)
                obj.env.Ygrid(:,:) = 0;
                obj.env.Ymin(:,:) = 0;
                obj.env.Ymax(:,:) = 0;
            elseif ((obj.env.tretire > 0) && (obj.env.tretire <= obj.env.T))
                obj.env.Ygrid(obj.env.tretire:obj.env.T,:) = 0;
                obj.env.Ymin(obj.env.tretire:obj.env.T,:) = 0;
                obj.env.Ymax(obj.env.tretire:obj.env.T,:) = 0;
            end
        end
        
        % Function to populate asset grid
        function AgridInit(obj)
            
            % Maximum assets
            obj.env.Agrid(1, obj.env.numPointsA) = obj.env.startA;
            for ixt = 2:obj.env.T+1
                obj.env.Agrid(ixt, obj.env.numPointsA) = (obj.env.Agrid(ixt-1, obj.env.numPointsA) + obj.env.Ymax(ixt-1,1) - obj.env.minCons)*(1 + obj.env.r);
            end
            
            % Minimum assets (including imposing borrowing constraint)
            obj.env.Agrid(obj.env.T+1, 1) = 0;
            for ixt = obj.env.T:-1:1
                obj.env.Agrid(ixt, 1) = obj.env.Agrid(ixt+1, 1)/(1+obj.env.r) - obj.env.Ymin(ixt,1) + obj.env.minCons;
                if ((obj.env.borrowingAllowed == 0) && (obj.env.Agrid(ixt, 1) < 0))
                    obj.env.Agrid(ixt, 1) = 0;
                end
            end
            
            % Need to check that maximum assets > minimum assets everywhere
            
            % Asset points in between
            for ixt = 1:obj.env.T+1
                %Agrid(ixt, :) = logspace(Agrid(ixt, 1), Agrid(ixt, obj.env.numPointsA), obj.env.numPointsA);
                % Next lines implement 3 log steps for comparison with
                % Cormac's code
                span = obj.env.Agrid(ixt, obj.env.numPointsA) - obj.env.Agrid(ixt, 1);
                loggrid = linspace(log(1+log(1+log(1))), log(1+log(1+log(1+span))), obj.env.numPointsA);
                grid = exp(exp(exp(loggrid)-1)-1)-1;
%                 % Next lines implement 5 log steps for comparison with
%                 % Cormac's code
%                 span = obj.env.Agrid(ixt, obj.env.numPointsA) - obj.env.Agrid(ixt, 1);
%                 loggrid = linspace(log(1+log(1+log(1+log(1+log(1))))), log(1+log(1+log(1+log(1+log(1+span))))), obj.env.numPointsA);
%                 grid = exp(exp(exp(exp(exp(loggrid)-1)-1)-1)-1)-1;
                obj.env.Agrid(ixt, :) = grid + obj.env.Agrid(ixt, 1)*ones(1, obj.env.numPointsA);

            end
        end
        
        % main solution function
        function solve(obj)
            if (obj.env.solveUsingEuler)
                solveUsingEuler(obj)
            else
                solveUsingValue(obj)
            end
        end
        
        % Solution maximising the value function rather than solving Euler
        % equation
        function solveUsingValue(obj)
            
            % Set terminal value function to zero
            obj.sln.value(obj.env.T+1, 1:obj.env.numPointsA) = 0;
            
            % Loop backwards in time
            for ixt = obj.env.T:-1:1
                
                % Loop across asset points
                for ixA = 1:1:obj.env.numPointsA
                    
                    % Information for optimisation
                    A0 = obj.env.Agrid(ixt, ixA);
                    Y = obj.env.Ygrid(ixt, 1);
                    lbA1 = obj.env.Agrid(ixt + 1, 1);
                    ubA1 = (A0 + Y - obj.env.minCons)*(1 + obj.env.r);
                    
                    % Need to check that ubA1 > lbA1
                    if (ubA1 - lbA1 < obj.env.tol)
                        negV = objectiveFunc(obj, ixt, A0, Y, lbA1);
                        A1 = lbA1;
                    else
                        % Compute solution
                        [A1, negV] = fminbnd(@(A1) objectiveFunc(obj, ixt, A0, Y, A1), lbA1, ubA1, optimset('TolX', obj.env.tol));
                    end
                        
                    % Store solution
                    obj.sln.value(ixt, ixA) = -negV;
                    obj.sln.policyA1(ixt, ixA) = A1;
                    obj.sln.policyC(ixt, ixA) = A0 + Y - A1/(1 + obj.env.r);
                    obj.sln.margUtil(ixt, ixA) = marginalUtility(obj, obj.sln.policyC(ixt, ixA));
                    
                end % ixA
                
            end % ixt
        end
        
        % Solution function based on Euler equation
        function solveUsingEuler(obj)

            % Set terminal (T+1) value function to zero
            obj.sln.value(obj.env.T+1, 1:obj.env.numPointsA) = 0;
            
            % Solve problem at T
            obj.sln.policyC(obj.env.T, :) = obj.env.Agrid(obj.env.T, :) + obj.env.Ygrid(obj.env.T, 1);
            obj.sln.policyA1(obj.env.T, :) = zeros(1, obj.env.numPointsA);
            obj.sln.value(obj.env.T, :) = utility(obj, obj.sln.policyC(obj.env.T, :));
            obj.sln.margUtil(obj.env.T, :) = marginalUtility(obj, obj.sln.policyC(obj.env.T, :));
            
            % Solve problem at T-1 back to 1
            for ixt = (obj.env.T-1):-1:1
                
                % Loop across asset grid points
                for ixA = 1:1:obj.env.numPointsA
                    
                    % Information for optimisation
                    A0 = obj.env.Agrid(ixt, ixA);
                    Y = obj.env.Ygrid(ixt, 1);
                    lbA1 = obj.env.Agrid(ixt + 1, 1);
                    ubA1 = (A0 + Y - obj.env.minCons)*(1 + obj.env.r);
                    
                    % Compute solution
                    lbSign = sign(eulerDifference(obj, ixt, A0, Y, lbA1));
                    % If liquidity constrained
                    if ((lbSign == 1) || (ubA1 - lbA1 < obj.env.tol))
                        obj.sln.policyA1(ixt, ixA) = lbA1;
                    else
                        ubSign = sign(eulerDifference(obj, ixt, A0, Y, ubA1));
                        if (lbSign*ubSign == 1)
                            error('Sign of Euler difference at lower bound and upper bound are the same. No solution to Euler equation');
                        end
                        obj.sln.policyA1(ixt, ixA) = fzero(@(A1) eulerDifference(obj, ixt, A0, Y, A1), [lbA1, ubA1], optimset('TolX', obj.env.tol));
                    end
                    
                    % Store solution
                    obj.sln.policyC(ixt, ixA) = A0 + Y - (obj.sln.policyA1(ixt, ixA) / (1 + obj.env.r));
                    obj.sln.value(ixt, ixA) = -objectiveFunc(obj, ixt, A0, Y, obj.sln.policyA1(ixt, ixA));
                    obj.sln.margUtil(ixt, ixA) = marginalUtility(obj, obj.sln.policyC(ixt, ixA));
                    
                end
                
            end
        
        end
        
        % Analytical policy functions
        function solveExact(obj)
            
            alpha = (obj.env.beta^(1/obj.env.gamma))*((1+obj.env.r)^((1-obj.env.gamma)/obj.env.gamma));
            
            for ixt = 1:1:obj.env.T
                periodsLeft = obj.env.T - ixt + 1;
                
                % Calculate discounted value of current and future income
                PVY = 0;
                discountFac = 1;
                for ixs= 0:1:(obj.env.T-ixt)
                    PVY = PVY + obj.env.Ygrid(ixt+ixs, 1) * discountFac;
                    discountFac = discountFac/(1 + obj.env.r);
                end
                
                for ixA = 1:1:obj.env.numPointsA
                    
                    % Calculate exact consumption policy function
                    if (abs(alpha - 1) < 1e-5)
                        obj.slnExact.policyC(ixt, ixA) = (obj.env.Agrid(ixt, ixA) + PVY) / periodsLeft;
                    else
                        obj.slnExact.policyC(ixt, ixA) = ((1 - alpha) / (1 - (alpha^periodsLeft))) * (obj.env.Agrid(ixt, ixA) + PVY);
                    end
                    
                    % Calculate exact value
                    discountFac = 1;
                    futureCons = obj.slnExact.policyC(ixt, ixA);
                    obj.slnExact.value(ixt, ixA) = utility(obj, futureCons);
                    for ixs = 1:1:obj.env.T-ixt
                        discountFac = discountFac * obj.env.beta;
                        futureCons = futureCons * ((obj.env.beta * (1 + obj.env.r)) ^ (1/obj.env.gamma));
                        obj.slnExact.value(ixt, ixA) = obj.slnExact.value(ixt, ixA) + (discountFac*utility(obj, futureCons));
                    end
                    
                    % Calculate exact marginal utility
                    obj.slnExact.margUtil(ixt, ixA) = marginalUtility(obj, obj.slnExact.policyC(ixt, ixA));
                    
                    % Calculate exact next-period assets policy function
                    obj.slnExact.policyA1(ixt, ixA) = (1 + obj.env.r) * (obj.env.Agrid(ixt, ixA) + obj.env.Ygrid(ixt, 1) - obj.slnExact.policyC(ixt, ixA));
                    
                end
            end
        end
        
        % Main simulation routine
        % Only need to simulate one individual because all individuals will
        % be identical
        function simulate(obj)
            
            obj.sims.a(1, 1) = obj.env.startA;
            
            for t = 1:1:obj.env.T
                obj.sims.v(t, 1) = interp1(obj.env.Agrid(t, :), obj.sln.value(t, :), obj.sims.a(t, 1), obj.env.interpMethod, 'extrap');
                obj.sims.a(t+1, 1) = interp1(obj.env.Agrid(t, :), obj.sln.policyA1(t, :), obj.sims.a(t, 1), obj.env.interpMethod, 'extrap');
                obj.sims.c(t, 1) = obj.sims.a(t, 1) + obj.env.Ygrid(t, 1) - (obj.sims.a(t+1, 1) / (1 + obj.env.r));
                obj.sims.du(t, 1) = marginalUtility(obj, obj.sims.c(t, 1));
            end
            
        end
        
        % Perform exact simulations
        function simulateExact(obj)
            
            alpha = (obj.env.beta^(1/obj.env.gamma))*((1+obj.env.r)^((1-obj.env.gamma)/obj.env.gamma));
            obj.simsExact.a(1, 1) = obj.env.startA;
            
            % Calculate PV of income
            PVY = 0;
            discountFac = 1;
            for t= 1:1:(obj.env.T)
                PVY = PVY + obj.env.Ygrid(t, 1) * discountFac;
                discountFac = discountFac/(1 + obj.env.r);
            end
            
            % Simulate consumption and assets
            for t = 1:1:obj.env.T
                obj.simsExact.c(t, 1) = (obj.simsExact.a(t, 1) + PVY) * (1 - alpha) / (1 - alpha^(obj.env.T - t + 1));
                obj.simsExact.du(t, 1) = marginalUtility(obj, obj.simsExact.c(t, 1));
                obj.simsExact.a(t+1, 1) = (obj.simsExact.a(t, 1) + obj.env.Ygrid(t, 1) - obj.simsExact.c(t, 1)) * (1 + obj.env.r);
                PVY = (PVY - obj.env.Ygrid(t, 1)) * (1 + obj.env.r);
            end
            
            % Simulate value
            obj.simsExact.v(obj.env.T, 1) = utility(obj, obj.simsExact.c(obj.env.T, 1));
            for t = obj.env.T-1:-1:1
                obj.simsExact.v(t, 1) = utility(obj, obj.simsExact.c(t, 1)) + (obj.env.beta * obj.simsExact.v(t+1, 1));
            end
        
        end
        
        
        % Plot analytical and numerical consumption policy functions
        function plotPolicyC(obj, year)
            figure(1)
            plot(obj.env.Agrid(year, :), obj.slnExact.policyC(year, :), 'g', 'LineWidth', 2)
            hold on;
            plot(obj.env.Agrid(year, :), obj.sln.policyC(year, :), 'r', 'LineWidth', 2)
            hold on;
            xlabel('Asset');
            ylabel('Policy (consumption)');
            legend('Analytical','Numerical', 4);
            title('Comparing numerical and analytical consumption functions')   
        end
        
        % Plot analytical and numerical next-period-asset policy functions
        function plotPolicyA1(obj, year)
            figure(2)
            plot(obj.env.Agrid(year, :), obj.slnExact.policyA1(year, :), 'g', 'LineWidth', 2)
            hold on;
            plot(obj.env.Agrid(year, :), obj.sln.policyA1(year, :), 'r', 'LineWidth', 2)
            hold on;
            xlabel('Asset');
            ylabel('Policy (next-period assets)');
            legend('Analytical','Numerical', 4);
            title('Comparing numerical and analytical next-period assets functions')   
        end

        % Plot analytical and numerical value functions
        function plotValue(obj, year)
            figure(3)
            plot(obj.env.Agrid(year, :), obj.slnExact.value(year, :), 'g', 'LineWidth', 2)
            hold on;
            plot(obj.env.Agrid(year, :), obj.sln.value(year, :), 'r', 'LineWidth', 2)
            hold on;
            xlabel('Asset');
            ylabel('Value function');
            legend('Analytical','Numerical', 4);
            title('Comparing numerical and analytical value functions')   
        end

        % Plot analytical and numerical marginal utility functions
        function plotMargUtil(obj, year)
            figure(4)
            plot(obj.env.Agrid(year, :), obj.slnExact.margUtil(year, :), 'g', 'LineWidth', 2)
            hold on;
            plot(obj.env.Agrid(year, :), obj.sln.margUtil(year, :), 'r', 'LineWidth', 2)
            hold on;
            xlabel('Asset');
            ylabel('Marginal utility function');
            legend('Analytical','Numerical', 4);
            title('Comparing numerical and analytical marginal utility functions')   
        end

        % Plot consumption
        function plotCons(obj)
            figure(5)
            plot(obj.simsExact.c, 'g', 'linewidth',2)
            hold on;
            plot(obj.sims.c, 'r', 'linewidth',2)
            hold on;
            xlabel('Age')
            ylabel('Consumption')
            legend('Analytical','Numerical', 4);
            title('Time path of consumption')
        end
        
        % Plot assets
        function plotAssets(obj)
            figure(6)
            plot(obj.simsExact.a, 'g', 'linewidth',2)
            hold on;
            plot(obj.sims.a, 'r', 'linewidth',2)
            hold on;
            xlabel('Age')
            ylabel('Assets')
            legend('Analytical','Numerical', 4);
            title('Time path of assets')
        end
        
        % Reset all the class properties to their initial empty state
        function clear(obj)
            obj.env = [];
            obj.sln = [];
            obj.slnExact = [];
            obj.sims = [];
            obj.simsExact = [];
        end
        
        
    end
end