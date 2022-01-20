classdef hayek
    %hayek Class with Hayek, 2016 model.
    %   With this class Hayek,2016 model of a vertical 1D infiltration on a
    %   soil is created. th can be plot and all variables get in h in time. 
    
    properties
        n= 3.5  %Exponent for relative permeability (kr=Se^n)
        a= 0.1  %Factor for WRC Se=Exp(a·h/n) [L-1]
        ksat = 1 %Saturated hydraulic conductivity [L·T-1]
        thres = 0.06 % Residual water content [-]
        thsat = 0.4  % Saturated water content [-]
        Sinit = 1.0  % Parameter to specify different initial conditions [-]
    end
    
    methods
        function obj = hayek(a,n,ksat,thres,thsat)
            %hayek Construct an instance of this class
            %   Detailed explanation goes here
            obj.a = a;
            obj.n = n;
            obj.thres = thres;
            obj.thsat = thsat;
            obj.ksat = ksat;
        end
        
        function outputArg = Sr(obj)
            %Sr Residual saturation.
            outputArg = obj.thres./obj.thsat;
        end

        function outputArg = kr_Se(obj,Se)
            %Sr Residual saturation.
            outputArg = Se.^(obj.n);
        end
        
        function outputArg = Se_h(obj,h)
            %Sr Residual saturation.
            outputArg = exp(obj.a.*h./obj.n);
        end      
        
        function outputArg = kr_h(obj,h)
            %Sr Residual saturation.
            outputArg =obj.kr_Se(obj.Se_h(h));
        end          
        
        function outputArg = h_Se(obj,Se)
            %Sr Residual saturation.
            outputArg = log(Se).*obj.n./obj.a;
        end            
        
        function outputArg = th_Se(obj,Se)
            %Sr Residual saturation.
            outputArg = obj.thres+(obj.thsat-obj.thres).*Se;
        end
        
        function outputArg = th_h(obj,h)
            %Sr Residual saturation.
            outputArg = obj.thres+(obj.thsat-obj.thres).*obj.Se_h(h);
        end
        
        %% Auxiliary expressions:
        
        function outputArg = K0(obj)
            %Sr Residual saturation.
            outputArg = obj.ksat/(obj.thsat-obj.thres);
        end
        
        function outputArg = v(obj)
            %Sr Residual saturation.
            outputArg = obj.K0().*obj.Sinit^(obj.n-1.0);
        end        
         
        function outputArg = chi0_z(obj,z)
            %Sr Residual saturation.
            outputArg = obj.chi_z_t(z,0.0);
        end   
        
        function outputArg = chi_z_t(obj,z,t)
            %Sr Residual saturation.
            outputArg = -z - obj.v*t;
        end           
 
        function outputArg = zf_t(obj,t)
            %Sr Residual saturation.
            outputArg = obj.chi0_z(0.0)+obj.v*t;
        end         
        
        function outputArg = Se_zvec_tsca(obj,z,t)
            %Sr Residual saturation.
            outputArg = zeros(size(z));
            zf_t = -obj.zf_t(t);
            outputArg(z>zf_t) = (1-exp(obj.a.*(obj.n-1).*(-z(z>zf_t)+zf_t)./obj.n)).^(1./(obj.n-1));
        end    
        
        function outputArg = h_zvec_tsca(obj,z,t)
            %Sr Residual saturation.
            outputArg = ln(obj.Se_zvec_tsca(z,t)).*obj.n./obj.a;           
        end
        
        function outputArg = plot_th_in_z_at_t(obj,z,t)
            %Sr Residual saturation.
            outputArg = plot(obj.th_Se(obj.Se_zvec_tsca(z,t)),z,'k');         
        end     
        
    end
end

