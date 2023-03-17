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

        measurement_uncertainty = 0.09; % FIXME calibration

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

        function v = get_acc_bias(obj)
            v = obj.x(7:9);
        end

        function P = get_process_cov(obj)
            P = obj.P;
        end


        function F = predict_jacobian(obj,vel_delta,dt)
            %{
            Returns the F matrix (n by n, n = len(x)) for EKF calculations.
            
            F = jacobian w/ respect to x of function f s.t.
            x^{.} = f(x^{^}, u, w)
            
            Note that the content of this matrix has been pre-determined,
            Only numerical values are set at run-time.
            %}

            dvn = vel_delta(1);
            dve = vel_delta(2);
            dvd = vel_delta(3);
            acc_bias_n = obj.x(7);
            acc_bias_e = obj.x(7);
            acc_bias_d = obj.x(7);
            
            F = [...
                1, 0, 0,    dt, 0, 0,                   0, 0, 0
                0, 1, 0,    0, dt, 0,                   0, 0, 0
                0, 0, 1,    0, 0, dt,                   0, 0, 0
                0, 0, 0,    (dvn - acc_bias_n), 0, 0,   1, 0, 0
                0, 0, 0,    0, (dve - acc_bias_e), 0,   0, 1, 0
                0, 0, 0,    0, 0, (dvd - acc_bias_d),   0, 0, 1
                0, 0, 0,    0, 0, 0,                    0, 0, 0
                0, 0, 0,    0, 0, 0,                    0, 0, 0
                0, 0, 0,    0, 0, 0,                    0, 0, 0];
        end

        function G = predict_process_noise(obj,w)
            %{
            Returns the G matrix (n by n, n = len(x)) for EKF calculations.
            
            G = jacobian w/ respect to w of function f s.t.
            x^{.} = f(x^{^}, u, w)
            
            Note that the content of this matrix has been pre-determined,
            Only numerical values are set at run-time.
            %}
            dvnCov = w(1);
            dveCov = w(2);
            dvdCov = w(3);

            G = [...
                0, 0, 0,    0, 0, 0,        0, 0, 0
                0, 0, 0,    0, 0, 0,        0, 0, 0
                0, 0, 0,    0, 0, 0,        0, 0, 0
                0, 0, 0,    dvnCov, 0, 0,   0, 0, 0
                0, 0, 0,    0, dveCov, 0,   0, 0, 0
                0, 0, 0,    0, 0, dvdCov,   0, 0, 0
                0, 0, 0,    0, 0, 0,        0, 0, 0
                0, 0, 0,    0, 0, 0,        0, 0, 0
                0, 0, 0,    0, 0, 0,        0, 0, 0];
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
            acc_bias_n = x(7);
            acc_bias_e = x(8);
            acc_bias_d = x(9);

            dvn = vel_delta(1);
            dve = vel_delta(2);
            dvd = vel_delta(3);

            x = [
                pn + vn * dt
                pe + ve * dt
                pd + vd * dt
                vn + (dvn - acc_bias_n)
                ve + (dve - acc_bias_e)
                vd + (dvd - acc_bias_d)
                acc_bias_n
                acc_bias_e
                acc_bias_d];
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
   
        function [x_new,P_new] = update_step_sensors(obj, z)
            %{
                "Meta-function", use this to call update
                Kept separate to split math and sensors.
            %}
            R = eye(3) * obj.measurement_uncertainty;
            h_x = sensors_measurement_model();
            H = sensors_measurement_jacobian();
            
            [x_new, P_new] = update_step(z, h_x, H, R);
        end

        function H = sensors_measurement_jacobian(obj)
            %{
                Jacobian matrix of the sensor model
            %}
            H = [...
                0, 0, 0,    1, 0, 0,    0, 0, 0
                0, 0, 0,    0, 1, 0,    0, 0, 0
                0, 0, 0,    0, 0, 1,    0, 0, 0];
        end

        function h_x = sensors_measurement_model(obj)
            %{
                Sensor model, currently : only uses velocity
            %}
            vn = obj.x(4);
            ve = obj.x(5);
            vd = obj.x(6);
            h_x = [vn ve vd]';
        end
       
        function obj = EKF_rocket(x_init,init_process_cov)
            obj.x = x_init;
            obj.P = ones(9)*init_process_cov;
        end

        function [obj,Qs,w] = set_additive_noise(obj,Ts)
            Fs = 1/Ts;

            obj.scale_var = 0.5*(1./(Fs.^2));
            obj.vel_delta_bias_sigma = obj.scale_var .* obj.AccelerometerBiasNoise;
            w = obj.scale_var.*[obj.AccelerometerNoise*ones(1,3)];

            Qs = diag([zeros(1, 3), obj.vel_delta_bias_sigma*ones(1,3), zeros(1,3)]);

        end


    end

end