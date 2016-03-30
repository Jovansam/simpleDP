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
            if (obj.env.gamma == 1)
                utils = log(cons);
            else
                utils = cons^(1 - obj.env.gamma) ./ (1 - obj.env.gamma);
            end
        end
        
        % Objective function
        % Rather than passing ixt (time index) as a separate argument, we
        % could pass a state structure combining ixt and A0
        function value = objectiveFunc(obj, ixt, A0, A1)
            cons = A0 - (A1 ./ (1 + obj.env.r));
            VA1 = interp1(obj.env.Agrid(ixt+1, :), obj.sln.value(ixt+1, :), A1, obj.env.interpMethod, 'extrap');
            value = utility(obj, cons) + (obj.env.beta * VA1);
            value = -value;
        end
    
    end
    methods (Access = public)
        
        % Function to initialise the environment
        function envInit(obj)
            obj.env.tol = 1e-10;
            obj.env.minCons = 1e-5;
            obj.env.T = 80;
            obj.env.r = 0.01;
            obj.env.beta = 1/(1+obj.env.r);
            obj.env.gamma = 1.5;
            obj.env.startA = 1;
            obj.env.numPointsA = 100;
            obj.env.interpMethod = 'pchip';  % Previously I had 'linear' but this gave very misleading simulations
            obj.env.Agrid = obj.AgridInit();
        end
        
        % Function to populate asset grid
        function Agrid = AgridInit(obj)
            
            % Maximum assets
            Agrid(1, obj.env.numPointsA) = obj.env.startA;
            for ixt = 2:obj.env.T+1
                Agrid(ixt, obj.env.numPointsA) = (Agrid(ixt-1, obj.env.numPointsA) - obj.env.minCons)*(1 + obj.env.r);
            end
            
            % Minimum assets
            Agrid(obj.env.T+1, 1) = 0;
            for ixt = obj.env.T:-1:1
                Agrid(ixt, 1) = Agrid(ixt+1, 1)/(1+obj.env.r) + obj.env.minCons;
            end
            
            % Need to check that maximum assets > minimum assets everywhere
            
            % Asset points in between
            for ixt = 1:obj.env.T+1
                %Agrid(ixt, :) = logspace(Agrid(ixt, 1), Agrid(ixt, obj.env.numPointsA), obj.env.numPointsA);
                % Next lines implement 5 log steps for comparison with
                % Cormac's code
                span = Agrid(ixt, obj.env.numPointsA) - Agrid(ixt, 1);
                loggrid = linspace(log(1+log(1+log(1+log(1+log(1))))), log(1+log(1+log(1+log(1+log(1+span))))), obj.env.numPointsA);
                grid = exp(exp(exp(exp(exp(loggrid)-1)-1)-1)-1)-1;
                Agrid(ixt, :) = grid + Agrid(ixt, 1)*ones(1, obj.env.numPointsA);

            end
        end

        % main solution function
        function solve(obj)
            
            % Set terminal value function to zero
            obj.sln.value(obj.env.T+1, 1:obj.env.numPointsA) = 0;
            
            % Loop backwards in time
            for ixt = obj.env.T:-1:1
                
                % Loop across asset points
                for ixA = 1:1:obj.env.numPointsA
                    
                    % Information for optimisation
                    A0 = obj.env.Agrid(ixt, ixA);
                    lbA1 = obj.env.Agrid(ixt + 1, 1);
                    ubA1 = (A0 - obj.env.minCons)*(1 + obj.env.r);
                    
                    % Need to check that ubA1 > lbA1
                    if (ubA1 - lbA1 < obj.env.tol)
                        negV = objectiveFunc(obj, ixt, A0, lbA1);
                        A1 = lbA1;
                    else
                        % Compute solution
                        [A1, negV] = fminbnd(@(A1) objectiveFunc(obj, ixt, A0, A1), lbA1, ubA1, optimset('TolX', obj.env.tol));
                    end
                        
                    % Store solution
                    obj.sln.value(ixt, ixA) = -negV;
                    obj.sln.policyA1(ixt, ixA) = A1;
                    obj.sln.policyC(ixt, ixA) = A0 - A1/(1 + obj.env.r);
                    
                end % ixA
                
            end % ixt
        end
       
        % Analytical policy functions
        % This is just copying Cormac's code
        function solveExact(obj)
            
            alpha = (obj.env.beta^(1/obj.env.gamma))*((1+obj.env.r)^((1-obj.env.gamma)/obj.env.gamma));
            
            for ixt = 1:1:obj.env.T
                periodsLeft = obj.env.T - ixt + 1;
                
                for ixA = 1:1:obj.env.numPointsA
                    
                    % Calculate exact consumption policy function
                    if (abs(alpha - 1) < 1e-5)
                        obj.slnExact.policyC(ixt, ixA) = obj.env.Agrid(ixt, ixA) / periodsLeft;
                    else
                        obj.slnExact.policyC(ixt, ixA) = ((1 - alpha) / (1 - (alpha^periodsLeft))) * obj.env.Agrid(ixt, ixA);
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
                    
                    % Calculate exact next-period assets policy function
                    obj.slnExact.policyA1(ixt, ixA) = (1 + obj.env.r) * (obj.env.Agrid(ixt, ixA) - obj.slnExact.policyC(ixt, ixA));
                    
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
                obj.sims.c(t, 1) = obj.sims.a(t, 1) - (obj.sims.a(t+1, 1) / (1 + obj.env.r));
            end
            
        end
        
        % Perform exact simulations
        function simulateExact(obj)
            
            alpha = (obj.env.beta^(1/obj.env.gamma))*((1+obj.env.r)^((1-obj.env.gamma)/obj.env.gamma));
            obj.simsExact.a(1, 1) = obj.env.startA;
            
            for t = 1:1:obj.env.T
                obj.simsExact.c(t, 1) = obj.simsExact.a(t, 1) * (1 - alpha) / (1 - alpha^(obj.env.T - t + 1));
                obj.simsExact.a(t+1, 1) = (obj.simsExact.a(t, 1) - obj.simsExact.c(t, 1)) * (1 + obj.env.r);
            end
            obj.simsExact.v(obj.env.T, 1) = utility(obj, obj.simsExact.c(obj.env.T, 1));
            for t = obj.env.T-1:-1:1
                obj.simsExact.v(t, 1) = utility(obj, obj.simsExact.c(t, 1)) + (obj.env.beta * obj.simsExact.v(t+1, 1));
            end
        
        end
        
        
        % Plot analytical and numerical consumption policy functions
        function plotPolicyC(obj, year)
            figure(1)
            plot(obj.env.Agrid(year, :), obj.sln.policyC(year, :), 'r', 'LineWidth', 2)
            hold on;
            plot(obj.env.Agrid(year, :), obj.slnExact.policyC(year, :), 'g', 'LineWidth', 2)
            hold on;
            xlabel('Asset');
            ylabel('Policy (consumption)');
            legend('Numerical','Analytical', 4);
            title('Comparing numerical and analytical consumption functions')   
        end
        
        % Plot analytical and numerical next-period-asset policy functions
        function plotPolicyA1(obj, year)
            figure(2)
            plot(obj.env.Agrid(year, :), obj.sln.policyA1(year, :), 'r', 'LineWidth', 2)
            hold on;
            plot(obj.env.Agrid(year, :), obj.slnExact.policyA1(year, :), 'g', 'LineWidth', 2)
            hold on;
            xlabel('Asset');
            ylabel('Policy (next-period assets)');
            legend('Numerical','Analytical', 4);
            title('Comparing numerical and analytical next-period assets functions')   
        end

        % Plot analytical and numerical value functions
        function plotValue(obj, year)
            figure(3)
            plot(obj.env.Agrid(year, :), obj.sln.value(year, :), 'r', 'LineWidth', 2)
            hold on;
            plot(obj.env.Agrid(year, :), obj.slnExact.value(year, :), 'g', 'LineWidth', 2)
            hold on;
            xlabel('Asset');
            ylabel('Value function');
            legend('Numerical','Analytical', 4);
            title('Comparing numerical and analytical value functions')   
        end
        
        % Plot consumption
        function plotCons(obj)
            figure(4)
            plot(obj.sims.c, 'r', 'linewidth',2)
            hold on;
            plot(obj.simsExact.c, 'g', 'linewidth',2)
            hold on;
            xlabel('Age')
            ylabel('Consumption')
            legend('Numerical','Analytical', 4);
            title('Time path of consumption')
        end
        
        % Plot assets
        function plotAssets(obj)
            figure(5)
            plot(obj.sims.a, 'r', 'linewidth',2)
            hold on;
            plot(obj.simsExact.a, 'g', 'linewidth',2)
            hold on;
            xlabel('Age')
            ylabel('Assets')
            legend('Numerical','Analytical', 4);
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