classdef Model < handle
    properties (Access = public)
        env
        sln
        sims
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
        function value = objectiveFunc(obj, ixt, ixY, A0, A1)
            Y = obj.env.Ygrid(ixt, ixY);
            cons = A0 + Y - (A1 ./ (1 + obj.env.r));
            VA1 = interp1(obj.env.Agrid(ixt+1, :), obj.sln.Evalue(ixt+1, :, ixY), A1, obj.env.interpMethod, 'extrap');
            value = utility(obj, cons) + (obj.env.beta * VA1);
            value = -value;
        end

        % Function for solving the Euler equation
        % Interpolates inverse marginal utility because this is linear
        function eulerDiff = eulerDifference(obj, ixt, ixY, A0, A1)
            Y = obj.env.Ygrid(ixt, ixY);
            invMargUtilAtA1 = interp1(obj.env.Agrid(ixt+1, :), inverseMarginalUtility(obj, obj.sln.EmargUtil(ixt+1, :, ixY)), A1, obj.env.interpMethod, 'extrap');
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
            obj.env.T = 40;
            obj.env.tretire = obj.env.T + 1;
            obj.env.r = 0.01;
            obj.env.beta = 0.98;
            obj.env.gamma = 1.5;
            obj.env.mu = 0;                     % Log income
            obj.env.rho = 0.75;
            obj.env.sigma = 0.25;
            obj.env.truncAt = 3;
            obj.env.startA = 0;
            obj.env.numPointsA = 20;
            obj.env.interpMethod = 'pchip';     % Previously I had 'linear' but this gave very misleading simulations
            obj.env.numPointsY = 8;
            obj.env.hcIncome = [0.1, 0.9]';
            obj.env.hcIncPDF = [0.5, 0.5]';
            obj.YgridInit();
            obj.env.borrowingAllowed = 1;
            obj.AgridInit();
            obj.env.solveUsingEuler = 0;
            obj.env.numSims = 2;
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
        function YgridInitOld(obj)
            obj.env.Ygrid = repmat(obj.env.hcIncome', obj.env.T, 1);
            obj.env.Ymin = min(obj.env.Ygrid, [], 2);
            obj.env.Ymax = max(obj.env.Ygrid, [], 2);
            obj.env.Q = repmat(obj.env.hcIncPDF', [obj.env.numPointsY, 1]);
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
        
        % Function to populate income grid
        function YgridInit(obj)
            
            % Process is lny = rho*lny_1 + (1-rho)*mu + epsilon
            % Mean is mu (ignores the truncation points)
            
            % Find std dev (ignores the truncation points)
            sigmalny = obj.env.sigma/((1-obj.env.rho^2)^0.5);
            
            % Find the N+1 boundary points and N expected values between
            % these boundary points

            % Initialise arrays
            lnY = NaN(obj.env.numPointsY, 1);
            lnYbounds = NaN(obj.env.numPointsY+1, 1);
            
            % Set lowest and highest bounds at truncation points
            lnYbounds(1,1) = obj.env.mu - (obj.env.truncAt*sigmalny);
            lnYbounds(obj.env.numPointsY+1, 1) = obj.env.mu + (obj.env.truncAt*sigmalny);
            
            % Find how much is excluded by truncation
            probAtLB = normcdf(lnYbounds(1, 1), obj.env.mu, sigmalny);
            probIn = normcdf(lnYbounds(obj.env.numPointsY+1, 1), obj.env.mu, sigmalny) - probAtLB;
            
            % Set the other bounds by splitting distribution into equal
            % probability areas
            for ixY = 2:1:obj.env.numPointsY
                lnYbounds(ixY,1) = norminv(probAtLB + probIn*(ixY-1)/obj.env.numPointsY, obj.env.mu, sigmalny);
            end

            % Find the N expected values between the boundary points
            % Calculate integral of x*pdf(x) between ixY_1 and ixY
            % fun is rescaled by numPointsY/probIn because expected value
            % is conditional on being in between lower and upper bounds
            fun = @(x) x .* normpdf(x,obj.env.mu,sigmalny) .* obj.env.numPointsY ./ probIn;
            for ixY = 1:1:obj.env.numPointsY
                lnY(ixY,1) = integral(fun, lnYbounds(ixY,1), lnYbounds(ixY+1,1));
            end

            % Find the probability of transitioning between different
            % ranges
            obj.env.Q = NaN(obj.env.numPointsY, obj.env.numPointsY);
            for ixi = 1:1:obj.env.numPointsY
                for ixj = 1:1:obj.env.numPointsY
                    minDraw = lnYbounds(ixj,1) - (1-obj.env.rho)*obj.env.mu - obj.env.rho*lnY(ixi,1);
                    maxDraw = lnYbounds(ixj+1,1) - (1-obj.env.rho)*obj.env.mu - obj.env.rho*lnY(ixi,1);
                    obj.env.Q(ixi,ixj) = normcdf(maxDraw/obj.env.sigma) - normcdf(minDraw/obj.env.sigma);
                end
                obj.env.Q(ixi,:) = obj.env.Q(ixi,:) ./ sum(obj.env.Q(ixi,:));
            end
            
            % Convert log income into income
            obj.env.Ygrid = repmat(exp(lnY'), obj.env.T, 1);
            obj.env.Ymin = repmat(exp(lnYbounds(1,1)), obj.env.T, 1);
            obj.env.Ymax = repmat(exp(lnYbounds(obj.env.numPointsY+1,1)), obj.env.T, 1);

            % Deal with retirement
            if (obj.env.tretire <= 0)
                obj.env.Ygrid(:,:) = 0;
                obj.env.Ymin(:,:) = 0;
                obj.env.Ymax(:,:) = 0;
            elseif ((obj.env.tretire > 0) && (obj.env.tretire <= obj.env.T))
                obj.env.Ygrid(obj.env.tretire:obj.env.T,:) = 0;
                obj.env.Ymin(obj.env.tretire:obj.env.T,:) = 0;
                obj.env.Ymax(obj.env.tretire:obj.env.T,:) = 0;
            end

            if (min(obj.env.Ygrid(:,1)) < 1e-4) || (max(obj.env.Ygrid(:,numPointsY)) > 1e5)
                warning('Combination of sigma and rho give a very high income variance. Numerical instability possible')
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
            
            % Set terminal value function and expected value function to zero
            obj.sln.value(obj.env.T+1, 1:obj.env.numPointsA, obj.env.numPointsY) = 0;
            obj.sln.Evalue(obj.env.T+1, 1:obj.env.numPointsA, obj.env.numPointsY) = 0;
            
            % Loop backwards in time
            for ixt = obj.env.T:-1:1
                
                % STEP 1: solve problem at grid points in assets and income
                
                % Loop across asset points
                for ixA = 1:1:obj.env.numPointsA
                    
                    % Loop across income grid points
                    for ixY = 1:1:obj.env.numPointsY
                        
                        % Information for optimisation
                        A0 = obj.env.Agrid(ixt, ixA);
                        Y = obj.env.Ygrid(ixt, ixY);
                        lbA1 = obj.env.Agrid(ixt + 1, 1);
                        ubA1 = (A0 + Y - obj.env.minCons)*(1 + obj.env.r);

                        % Need to check that ubA1 > lbA1
                        if (ubA1 - lbA1 < obj.env.tol)
                            negV = objectiveFunc(obj, ixt, ixY, A0, lbA1);
                            A1 = lbA1;
                        else
                            % Compute solution
                            [A1, negV] = fminbnd(@(A1) objectiveFunc(obj, ixt, ixY, A0, A1), lbA1, ubA1, optimset('TolX', obj.env.tol));
                        end

                        % Store solution
                        obj.sln.value(ixt, ixA, ixY) = -negV;
                        obj.sln.policyA1(ixt, ixA, ixY) = A1;
                        obj.sln.policyC(ixt, ixA, ixY) = A0 + Y - A1/(1 + obj.env.r);
                        obj.sln.margUtil(ixt, ixA, ixY) = marginalUtility(obj, obj.sln.policyC(ixt, ixA, ixY));
                    
                    end % ixY
                    
                    % STEP 2: integrate out income today conditional on
                    % income yesterday to get EV and EdU
                    
                    % Loop across income grid points
                    for ixY = 1:1:obj.env.numPointsY
                        % In these Evalue and EmargUtil matrices, ixA
                        % refers to time ixt while ixY refers to time ixt-1
                        obj.sln.Evalue(ixt, ixA, ixY) = obj.env.Q(ixY, :) * squeeze(obj.sln.value(ixt, ixA, :));
                        obj.sln.EmargUtil(ixt, ixA, ixY) = obj.env.Q(ixY, :) * squeeze(obj.sln.margUtil(ixt, ixA, :));
                    end % ixY
                    
                end % ixA
                
                fprintf('Period %d complete\n', ixt);
                
            end % ixt
        end
        
        % Solution function based on Euler equation
        function solveUsingEuler(obj)

            % Set terminal (T+1) value functions to zero
            obj.sln.value(obj.env.T+1, 1:obj.env.numPointsA, obj.env.numPointsY) = 0;
            obj.sln.Evalue(obj.env.T+1, 1:obj.env.numPointsA, obj.env.numPointsY) = 0;
            
            % Solve problem at T
            for ixA = 1:1:obj.env.numPointsA
                for ixY = 1:1:obj.env.numPointsY
                    obj.sln.policyC(obj.env.T, ixA, ixY) = obj.env.Agrid(obj.env.T, ixA) + obj.env.Ygrid(obj.env.T, ixY);
                    obj.sln.policyA1(obj.env.T, ixA, ixY) = 0;
                    obj.sln.value(obj.env.T, ixA, ixY) = utility(obj, obj.sln.policyC(obj.env.T, ixA, ixY));
                    obj.sln.margUtil(obj.env.T, ixA, ixY) = marginalUtility(obj, obj.sln.policyC(obj.env.T, ixA, ixY));
                end
                for ixY = 1:1:obj.env.numPointsY
                    obj.sln.Evalue(obj.env.T, ixA, ixY) = obj.env.Q(ixY, :) * squeeze(obj.sln.value(obj.env.T, ixA, :));
                    obj.sln.EmargUtil(obj.env.T, ixA, ixY) = obj.env.Q(ixY, :) * squeeze(obj.sln.margUtil(obj.env.T, ixA, :));
                end
            end
            fprintf('Period %d complete\n', obj.env.T);

            % Solve problem at T-1 back to 1
            for ixt = (obj.env.T-1):-1:1
                
                % Loop across asset grid points
                for ixA = 1:1:obj.env.numPointsA
                    
                    % Loop across income grid points
                    for ixY = 1:1:obj.env.numPointsY
                    
                        % Information for optimisation
                        A0 = obj.env.Agrid(ixt, ixA);
                        Y = obj.env.Ygrid(ixt, ixY);
                        lbA1 = obj.env.Agrid(ixt + 1, 1);
                        ubA1 = (A0 + Y - obj.env.minCons)*(1 + obj.env.r);

                        % Compute solution
                        lbSign = sign(eulerDifference(obj, ixt, ixY, A0, lbA1));
                        % If liquidity constrained
                        if ((lbSign == 1) || (ubA1 - lbA1 < obj.env.tol))
                            obj.sln.policyA1(ixt, ixA, ixY) = lbA1;
                        else
                            ubSign = sign(eulerDifference(obj, ixt, ixY, A0, ubA1));
                            if (lbSign*ubSign == 1)
                                error('Sign of Euler difference at lower bound and upper bound are the same. No solution to Euler equation');
                            end
                            obj.sln.policyA1(ixt, ixA, ixY) = fzero(@(A1) eulerDifference(obj, ixt, ixY, A0, A1), [lbA1, ubA1], optimset('TolX', obj.env.tol));
                        end

                        % Store solution
                        obj.sln.policyC(ixt, ixA, ixY) = A0 + Y - (obj.sln.policyA1(ixt, ixA, ixY) / (1 + obj.env.r));
                        obj.sln.value(ixt, ixA, ixY) = -objectiveFunc(obj, ixt, ixY, A0, obj.sln.policyA1(ixt, ixA, ixY));
                        obj.sln.margUtil(ixt, ixA, ixY) = marginalUtility(obj, obj.sln.policyC(ixt, ixA, ixY));
                        
                    end
                    
                    % Loop across income grid points
                    for ixY = 1:1:obj.env.numPointsY
                        obj.sln.Evalue(ixt, ixA, ixY) = obj.env.Q(ixY, :) * squeeze(obj.sln.value(ixt, ixA, :));
                        obj.sln.EmargUtil(ixt, ixA, ixY) = obj.env.Q(ixY, :) * squeeze(obj.sln.margUtil(ixt, ixA, :));
                    end
                    
                end
                fprintf('Period %d complete\n', ixt);

            end
        
        end
        
       
        % Simulate
        function simulate(obj)
            
            % Initialise arrays that will hold simulations
            obj.sims.y = NaN(obj.env.T, obj.env.numSims);
            obj.sims.c = NaN(obj.env.T, obj.env.numSims);
            obj.sims.v = NaN(obj.env.T, obj.env.numSims);
            obj.sims.a = NaN(obj.env.T+1, obj.env.numSims);
            ixY = NaN(obj.env.T, obj.env.numSims);
            
            % STEP 1: obtain path for exogenous income
            for t = 1:1:obj.env.T
                seed = t;
                [obj.sims.y(t, :), ixY(t, :)] = getDiscreteDraws(obj.env.Ygrid(t, :), ...
                    obj.env.hcIncPDF, 1, obj.env.numSims, seed);
                if (t >= obj.env.tretire)
                    ixY(t, :) = 1;
                end
            end
            
            % STEP 2: simulate endogenous outcomes
            for s = 1:1:obj.env.numSims
                obj.sims.a(1, s) = obj.env.startA;
                for t = 1:1:obj.env.T
                    obj.sims.a(t+1, s) = interp1(obj.env.Agrid(t, :), obj.sln.policyA1(t, :, ixY(t, s)), ...
                        obj.sims.a(t, s), obj.env.interpMethod, 'extrap');
                    obj.sims.c(t, s) = obj.sims.a(t, s) + obj.sims.y(t, s) - ...
                        (obj.sims.a(t+1, s) / (1 + obj.env.r));
                    obj.sims.v(t, s) = interp1(obj.env.Agrid(t, :), obj.sln.value(t, :, ixY(t, s)), ...
                        obj.sims.a(t, s), obj.env.interpMethod, 'extrap');
                end % t
            end % s
            
        end
        
        
        
        % Plot analytical and numerical consumption policy functions
        function plotPolicyC(obj, year)
            figure(1)
            plot(obj.env.Agrid(year, :), obj.sln.policyC(year, :, obj.env.numPointsY), 'g', 'LineWidth', 2)
            hold on;
            plot(obj.env.Agrid(year, :), obj.sln.policyC(year, :, 1), 'r', 'LineWidth', 2)
            hold on;
            xlabel('Asset');
            ylabel('Policy function (consumption function)');
            legend('Higest income', 'Lowest income', 4);
            title('Policy function (consumption function)')
        end

        % Plot analytical and numerical next-period-asset policy functions
        function plotPolicyA1(obj, year)
            figure(2)
            plot(obj.env.Agrid(year, :), obj.sln.policyA1(year, :, obj.env.numPointsY), 'g', 'LineWidth', 2)
            hold on;
            plot(obj.env.Agrid(year, :), obj.sln.policyA1(year, :, 1), 'r', 'LineWidth', 2)
            hold on;
            xlabel('Asset');
            ylabel('Policy function (next-period assets)');
            legend('Higest income', 'Lowest income', 4);
            title('Policy function (next-period assets)')
        end

        % Plot analytical and numerical value functions
        function plotValue(obj, year)
            figure(3)
            plot(obj.env.Agrid(year, :), obj.sln.value(year, :, obj.env.numPointsY), 'g', 'LineWidth', 2)
            hold on;
            plot(obj.env.Agrid(year, :), obj.sln.value(year, :, 1), 'r', 'LineWidth', 2)
            hold on;
            xlabel('Asset');
            ylabel('Value function');
            legend('Higest income', 'Lowest income', 4);
            title('Value function')   
        end
        
        % margUtil Evalue EmargUtil

        % Plot analytical and numerical marginal utility functions
        function plotMargUtil(obj, year)
            figure(4)
            plot(obj.env.Agrid(year, :), obj.sln.margUtil(year, :, obj.env.numPointsY), 'g', 'LineWidth', 2)
            hold on;
            plot(obj.env.Agrid(year, :), obj.sln.margUtil(year, :, 1), 'r', 'LineWidth', 2)
            hold on;
            xlabel('Asset');
            ylabel('Marginal utility function');
            legend('Higest income', 'Lowest income', 4);
            title('Marginal utility function')   
        end

        % Plot consumption
        function plotCons(obj)
            figure(5)
            plot(obj.sims.c, 'linewidth',2)
            xlabel('Age')
            ylabel('Consumption')
            title('Time path of consumption')
        end
        
        % Plot assets
        function plotAssets(obj)
            figure(6)
            plot(obj.sims.a, 'linewidth',2)
            xlabel('Age')
            ylabel('Assets')
            title('Time path of assets')
        end
        
        % Reset all the class properties to their initial empty state
        function clear(obj)
            obj.env = [];
            obj.sln = [];
            obj.sims = [];
        end
        
        
    end
end