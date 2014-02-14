classdef noise_and_error
    % This class contains the noise learning and error estimations
    
    properties
        true_error; convergence; var_noise_old; var_noise; conv_crit; dump_learn;
    end
    
    methods
        
        % Constructor function
        function obj = noise_and_error(conv_crit,var_noise_init,dump_learn)
            obj.conv_crit = conv_crit; obj.var_noise_old = var_noise_init; obj.var_noise = var_noise_init; obj.dump_learn = dump_learn;
        end
        
        % learning of the noise by expectation maximisation
        function obj = learn_noise(obj,Y,W,V)
            obj.var_noise = obj.dump_learn .* obj.var_noise_old + (1 - obj.dump_learn) .* ((Y - W).^2 * (1 + V ./ obj.var_noise_old).^(-2).' ) ./ sum((1 + V ./ obj.var_noise_old).^(-1) );
            if (obj.var_noise < 1e-100); obj.var_noise = 1e-100; end;
            obj.var_noise_old = obj.var_noise;
        end
        
        % true MSE (with respect to the original signal)
        function obj = compute_true_MSE(obj,signal,X)
            if (isreal(signal) ); obj.true_error = sum((signal - X).^2) ./ max(size(X) );
            else obj.true_error = 0.5 .* sum((real(signal) - real(X) ).^2 + (imag(signal) - imag(X) ).^2) ./ max(size(X) ); end
            
        end
        
        % section error rate for superposition codes
        function obj = compute_true_SER(obj,signal,X,varargin)
            if (max(max(X) ) > 0)
                Xtrue = reshape(X, varargin{1}, [] );
                [Y, I] = max(Xtrue);
                BB = zeros(size(Xtrue) );
                BB(sub2ind(size(Xtrue), I, 1 : length(I) ) ) = 1;
                obj.true_error = max(size(X) ) / varargin{1};
                BB = reshape(BB, size(X) );
            else
                BB = X;
            end
            obj.true_error = max(size(X) ) / varargin{1};
            for l = 1 : max(size(X) ) / varargin{1}
                obj.true_error = obj.true_error - isequal(signal((l - 1) * varargin{1} + 1 : l * varargin{1} ), varargin{2}(l) .* BB((l - 1) * varargin{1} + 1 : l * varargin{1} ) );
            end
            obj.true_error = obj.true_error ./ max(size(X) ) .* varargin{1};
        end
        
        % convergence
        function obj = compute_convergence(obj,X_old,X)
            obj.convergence = sum(vabs(X_old - X) ) ./ max(size(X) );
        end
        
    end
    
end