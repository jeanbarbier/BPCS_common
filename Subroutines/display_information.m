% shows the first line of information

if (opt.print > 0)
    switch opt.prior
        case 'SparseGauss'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   m  var   damping   FreeEntropy   convergence')
            else
                disp('iter   rho   noise   m  var   damping   FreeEntropy   convergence   error')
            end
        case 'SparseGaussPositive'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   m  var   damping   FreeEntropy   convergence')
            else
                disp('iter   rho   noise   m  var   damping   FreeEntropy   convergence   error')
            end
        case 'SparseGaussCut'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   m  var   damping   FreeEntropy   convergence')
            else
                disp('iter   rho   noise   m  var   damping   FreeEntropy   convergence   error')
            end
        case '2Gauss'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   m_1   m_2   var_1   var_2   damping   FreeEntropy   convergence')
            else
                disp('iter   rho   noise   m_1   m_2   var_1   var_2   damping   FreeEntropy   convergence   error')
            end
        case 'SparseExponential'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   expo   damping   FreeEntropy   convergence')
            else
                disp('iter   rho   noise   expo   damping   FreeEntropy   convergence   error')
            end
        case 'SparseConstant'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   c_down   c_up   damping   FreeEntropy   convergence')
            else
                disp('iter   rho   noise   c_down   c_up   damping   FreeEntropy   convergence   error')
            end
        case 'SparseBinary'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   damping   FreeEntropy   convergence')
            else
                disp('iter   rho   noise   damping   FreeEntropy   convergence   error')
            end
        case 'L1'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   damping   FreeEntropy   convergence')
            else
                disp('iter   rho   noise   damping   FreeEntropy   convergence   error')
            end
        case 'Laplace'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   beta   damping   FreeEntropy   convergence')
            else
                disp('iter   rho   noise   beta   damping   FreeEntropy   convergence   error')
            end
        case 'Binary1'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   damping   FreeEntropy   convergence')
            else
                disp('iter   rho   noise   damping   FreeEntropy   convergence   error')
            end
        case 'Complex'
            if (max(size(opt.signal) ) < 2)
                disp('iter   rho   noise   m  var   damping   FreeEntropy   convergence')
            else
                disp('iter   rho   noise   m  var   damping   FreeEntropy   convergence   error')
            end
        otherwise
            disp('unknown prior')
    end
end