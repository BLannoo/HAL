%% Script containing a collection of tools for GRN analysis
classdef GRN_tools
    methods(Static)
        function params = Cube_sampling(kInitial,numSamples,scale)
            % Dimention of the space to be sampled
            dim = length(kInitial);

            % random values between -1 and 1
            r = rand(dim,numSamples)*2-1;

            % params in a logaritmic cube aroud kInitial
            params = exp(r*log(scale)+log(kInitial)'*ones(1,numSamples))';
        end

        function [param,exponents,spearmanCorr] = PCA(k,per)
            if sum(per<=0) > 0
                warning(['All data should have a period bigger then 0,'...
                    ' did you remove non-oscillatory data points?'])
            end

            data = [k' ; per];

            N = length(per);

            % Moving analysis to the log-scale
            data = log(data);

            % Calculating the covariance matrix
            covarianceMat = cov(data');

            % Getting the smallest eigenvector
            [V,D] = eig(covarianceMat);
            exponents = D(1,1).*V(1:end-1,1);

            % Calculating the predicting parameter (previously known as beta)
            param = prod(k .^ (ones(N,1) * exponents'),2);

            % Evaluating the predicting strenght of the resulting parameter
            spearmanCorr = corr(param,per','type','Spearman');
        end

    end
end