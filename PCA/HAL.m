%% Script containing a collection of functions connected to the HAL model
classdef HAL<handle
    properties %(SetAccess = private)
        % Defining
        k % parameters
        
        % Settings
        x0 % initial state
        tEnd % length of the simulation
        
        % Flags
        toFewPer
        
        % Derived
        numPer
        phaseB % length of the B phase
        phaseA % length of the A phase
        maxA % max value of A (after transient)
        minA % min value of A (after transient)
        maxB % max value of A (after transient)
        minB % min value of A (after transient)
    end
    properties (SetAccess = private, Transient)
        track % solution of the set
    end
    properties (Dependent)
        beta
        period
    end
    methods
        %% A standard constructor with fixed k
        function obj = HAL(k)
            obj.k = k;
            % 1 gene and 1 A protein
            % (starting A at 1 garanties A phase comes first)
            obj.x0 = [1 0 1 0 0];
            % standard simulation length
            obj.tEnd = 100;
        end
        %% Evaluate HAL
        function eval(obj)
            % 0 periods known originally
            obj.toFewPer = 1;
            % Save the original tEnd to limit
            tEndOriginal = obj.tEnd;
            while obj.toFewPer && obj.tEnd < tEndOriginal*50
                obj.solve;
                obj.calcPeriod;
                if obj.toFewPer
                    obj.tEnd = obj.tEnd*2;
                end
            end
            
            obj.calcAmp;
        end
        %% Solve this equation
        function solve(obj)
            equ = @(x,k)[k(1)*(1-x(1)) - k(7)*x(1)*x(3) ; ...
                k(3)*x(1) + k(8) * (1-x(1)) - k(2)*x(2) ; ...
                k(4)*x(2) - k(9)*x(3) - k(6)*x(3)*x(4) + k(12)*x(5) ...
                    + k(1)*(1-x(1)) - k(7)*x(1)*x(3) ; ...
                k(5) - k(10) * x(4) - k(6)*x(3)*x(4) + k(12)*x(5) ; ...
                k(6)*x(3)*x(4) - k(12)*x(5) - k(11)*x(5)];
            if isempty(obj.track)
                obj.track = ode23tb(@(t,x)equ(x,obj.k),[0 obj.tEnd],obj.x0);
            else
                obj.track = odextend(obj.track,[],obj.tEnd);
            end
        end
        %% Get beta
        function beta = get.beta(obj)
            beta = obj.k(1)*obj.k(3)*obj.k(4)/(obj.k(5)*obj.k(2)^2);
        end
        %% Get period
        function period = get.period(obj)
            period = obj.phaseA + obj.phaseB;
        end 
    end
    methods(Access=private)
        %% Calculate period
        function calcPeriod(obj)
            % 0 => A phase; 1 => B phase
            phase = obj.track.y(3,:) < obj.track.y(4,:);
            % list of transition times from A phase to B phase
            AtoB = obj.track.x(diff(phase)==1);
            % list of transition times from B phase to A phase
            BtoA = obj.track.x(diff(phase)==-1);

            % if unequal periods of A and B remove one
            if length(BtoA) < length(AtoB)
                AtoB(end) = [];
            end

            obj.numPer = length(AtoB);
            
            % break the evaluation if to little periods found
            if(obj.numPer<10)
                obj.toFewPer = 1;
                return
            else
                obj.toFewPer = 0;
            end

            % estimate the length of the B phase by the last 5 instances
            Bphases = BtoA-AtoB;
            obj.phaseB = mean(Bphases(end-5:end));

            % estimate the length of the A phase by the last 5 instances
            Aphases = AtoB(2:end)-BtoA(1:end-1);
            obj.phaseA = mean(Aphases(end-5:end));
        end
        %% Calculate amplitude of the limit cycle
        function calcAmp(obj)
            if obj.numPer < 10
                warning(['calculate amplitude called on a track with'...
                    ' less then 10 periods. This might not give a'...
                    ' correct result.'])
            end
            
            % Calculate maxima and minima of the limit cycle for A and B
            % (last 4/5-th to skip the transient)
            obj.maxA = max(obj.track.y(3,floor(end*4/5):end));
            obj.minA = min(obj.track.y(3,floor(end*4/5):end));
            obj.maxB = max(obj.track.y(4,floor(end*4/5):end));
            obj.minB = min(obj.track.y(4,floor(end*4/5):end));
        end
    end
end