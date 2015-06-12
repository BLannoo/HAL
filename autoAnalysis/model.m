classdef model < handle
    properties
        k % the parameters of the model
        x0 % the initial state of the model
        
        % Derived quantities
        oscil % boolean whether oscillatory
        numPer % number of periods
        phaseA % length of the A phase
        phaseB % length of the B phase
        maxA % max value of A (after transient)
        minA % min value of A (after transient)
        maxB % max value of A (after transient)
        minB % min value of A (after transient)
    end
    properties(Dependent)
        beta % theta * uA / uB
        period % the total period of the limit cycle oscillations
    end
    methods
        %% Get beta
        function beta = get.beta(obj)
            beta = obj.k(1)*obj.k(3)*obj.k(4)/(obj.k(5)*obj.k(2)^2);
        end
        %% Get period
        function period = get.period(obj)
            period = obj.phaseA + obj.phaseB;
        end
    end
    methods
        % Set all the properties to standard values
        function obj = model()
            obj.k = zeros(1,12);
            obj.x0 = [1 0 1 0 0];
        end
    end
end