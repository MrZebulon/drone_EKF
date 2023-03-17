classdef EKF_rocket
    properties
        %{
            x: system state
            P: system covariance
            Ts: integration period
        %}
        
        x;
        P;
        Ts;

        AccelerometerNoise = 2;

        % extra additive noise
        AccelerometerBiasNoise =  2e-4;
        additiveNoise = 1e-8;

        scale_var = -1;
        vel_delta_bias_sigma = -1;

    end

    methods

        function obj = set_state(obj,x)
            obj.x = x;
        end

        function obj = set_process_cov(obj,P)
            obj.P = P;
        end

        function x = get_state(obj)
            x = obj.x;
        end

        function p = get_posNED(obj)
            p = obj.x(1:3);
        end

        function v = get_velNED(obj)
            v = obj.x(4:6);
        end

        function P = get_process_cov(obj)
            P = obj.P;
        end


        function F = predict_jacobian(obj,vel_delta,dt)
            %{
            Returns the F matrix (n by n, n = len(x)) for EKF calculations.
            
            F = jacobian w/ respect to x of function f s.t.
            x^{.} = f(x^{^}, u) + w
            
            Note that the content of this matrix has been pre-determined,
            Only numerical values are set at run-time.
            %}
            
            F = [...
                0, 0, 0, dt, 0, 0
                0, 0, 0, 0, dt, 0
                0, 0, 0, 0, 0, dt
                0, 0, 0, 0, 0, 0
                0, 0, 0, 0, 0, 0
                0, 0, 0, 0, 0, 0];
        end

        function G = predict_process_noise(obj,w)
            %{
            Returns the G matrix (n by n, n = len(x)) for EKF calculations.
            
            G = jacobian w/ respect to w of function f s.t.
            x^{.} = f(x^{^}, u) + w
            
            Note that the content of this matrix has been pre-determined,
            Only numerical values are set at run-time.
            %}

            G = [0, 0, 0, 0, 0, 0
                0, 0, 0, 0, 0, 0
                0, 0, 0, 0, 0, 0
                0, 0, 0, 0, 0, 0
                0, 0, 0, 0, 0, 0
                0, 0, 0, 0, 0, 0];
        end

        function x = predict_state(obj,vel_delta,dt)
            %{
                Realises prediction computations.
                Returns the a priori prediction. 
            %}
            x = obj.x;

            pn = x(1);
            pe = x(2);
            pd = x(3);
            vn = x(4);
            ve = x(5);
            vd = x(6);

            x = [
                pn + vn * dt
                pe + ve * dt
                pd + vd * dt
                vn
                ve
                vd];
        end

        function obj = predict_step(obj,accB,Ts)
            %{
                Realises all the steps for a priori estimation.
            %}
            
            % state prediction

            vel_delta = accB'*Ts;

            x_new = obj.predict_state(vel_delta,Ts);

            % covariance prediction

            [~,Qs,w] = obj.set_additive_noise(Ts);

            G = obj.predict_process_noise(w);
            F = obj.predict_jacobian(vel_delta,Ts);

            P_new = F*obj.P*(F')+G+Qs;
            P_new = 0.5*(P_new+P_new');

            obj.x = x_new;
            obj.P = P_new;

        end

        function [x_new,P_new]  = update_step(obj,z,h_x,H,R)
            %{
                Turns the a priori prediction into the a posteriori one for a given measurement point
            %}
            
            nx = size(obj.x,1);
            inov = z-h_x;
            S = H*obj.P*(H')+R;
            K = obj.P*(H')*inv(S);
            x_new = obj.x + K*inov;
            P_new = (eye(nx)-K*H)*obj.P;
        end
   
        function [x_new,P_new] = update_step_sensors(x, P, z)
        
            R = eye(3) * 0; % FIXME calibration
            h_x = sensors_measurement_model();
            H = sensors_measurement_jacobian();
            
            [x_new, P_new] = update_step(x, P, z, h_x, H, R);
        end

        function H = sensors_measurement_jacobian(obj)
            H = [...
                0, 0, 0, 1, 0, 0
                0, 0, 0, 0, 1, 0
                0, 0, 0, 0, 0, 1];
        end

        function h_x = sensors_measurement_model(obj)
            x = obj.x;
            vn = x(4);
            ve = x(5);
            vd = x(6);
            h_x = [vn ve vd]'; % FIXME ?
        end
       
        function obj = EKF_rocket(x_init,init_process_cov)
            obj.x = x_init;
            obj.P = ones(22)*init_process_cov;
        end

        function [obj,Qs,w] = set_additive_noise(obj,Ts)
            Fs = 1/Ts;

            obj.scale_var = 0.5*(1./(Fs.^2));
            obj.vel_delta_bias_sigma = obj.scale_var .* obj.AccelerometerBiasNoise;

            w = obj.scale_var.*[obj.AccelerometerNoise*ones(1,3)];

            Qs = diag([obj.additiveNoise.*ones(1,10), obj.vel_delta_bias_sigma*ones(1,3)]);

        end


    end

end