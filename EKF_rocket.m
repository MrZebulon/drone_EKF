classdef EKF_rocket
    %
    % State vector definition
    %   1-3 : position
    %   4-6 : speed
    %   7-9 : accel bias
    %   10 : baro bias
    %
    % Control vector definition
    %   1-3 : acceleration
    %
    % Noise vector definition
    %   1-3 : acceleration noise
    %   4 : barometer noise
    %   5-7 : acceleration bias noise
    %   8 : barometer bias noise
    properties
        x;
        P;
        Ts;

        % units = m/s
        accel_bias = [-0.026842656242568 0.033420780321046 -0.007947030636161];
        accel_noise = 1;
        accel_bias_noise = 2e-4;

        % units = m
        baro_bias = 361.3487972164834;
        baro_noise = 0.0011529662809109404;
        baro_bias_noise = 2.8486463440220755e-06;

        baro_measurement_uncertainty = 0.1;
        
        additiveNoise = 1e-8;

        scale_var = -1;
        vel_delta_bias_sigma = -1;
        pos_delta_bias_sigma = -1;


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

        function b_acc = get_acc_bias(obj)
            b_acc = obj.x(7:9);
        end

        function b_baro = get_baro_bias(obj)
            b_baro = obj.x(10);
        end

        function P = get_process_cov(obj)
            P = obj.P;
        end

        function F = predict_jacobian_preintegrated(obj, u, dt)
            %{
                Returns a preintegrated F matrix (jacobian w/ respect to x of
                function f)
            %}
            
             F = [...
                0, 0, 0,    obj.dt, 0, 0,    0, 0, 0,      0
                0, 0, 0,    0, obj.dt, 0,    0, 0, 0,      0
                0, 0, 0,    0, 0, obj.dt,    0, 0, 0,      0
                0, 0, 0,    0, 0, 0,         1, 0, 0,      0
                0, 0, 0,    0, 0, 0,         0, 1, 0,      0
                0, 0, 0,    0, 0, 0,         0, 0, 1,      0
                0, 0, 0,    0, 0, 0,         0, 0, 0,      0
                0, 0, 0,    0, 0, 0,         0, 0, 0,      0
                0, 0, 0,    0, 0, 0,         0, 0, 0,      0
                0, 0, 0,    0, 0, 0,         0, 0, 0,      0];

             F = F + eye(size(F));
        end

        function Q = predict_covariance_preintegrated(obj, w)
            % The Q matrix is the process Covariance Matrix.
            % Q(i,j) = Cov(x_i, x_j)

            dvxCov = w(1);
            dvyCov = w(2);
            dvzCov = w(3);
                       
           Q = [...
                0, 0, 0,    0, 0, 0,                    0, 0, 0,      0
                0, 0, 0,    0, 0, 0,                    0, 0, 0,      0
                0, 0, 0,    0, 0, 0,                    0, 0, 0,      0
                0, 0, 0,    dvxCov*obj.Ts^2, 0, 0,      0, 0, 0,      0
                0, 0, 0,    0, dvyCov*obj.Ts^2, 0,      0, 0, 0,      0
                0, 0, 0,    0, 0, dvzCov*obj.Ts^2,      0, 0, 0,      0
                0, 0, 0,    0, 0, 0,                    0, 0, 0,      0
                0, 0, 0,    0, 0, 0,                    0, 0, 0,      0
                0, 0, 0,    0, 0, 0,                    0, 0, 0,      0
                0, 0, 0,    0, 0, 0,                    0, 0, 0,      0];
        end

        function x = predict_state(obj,u,dt)
            %{
                Performs the a priori prediction.
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

            bias_baro = x(10);

            an = u(1);
            ae = u(2);
            ad = u(3);

            x = [
                pn + vn * dt
                pe + ve * dt
                pd + vd * dt
                vn + (an + acc_bias_n) * dt
                ve + (ae + acc_bias_e) * dt
                vd + (ad + acc_bias_d) * dt
                acc_bias_n
                acc_bias_e
                acc_bias_d
                bias_baro];
        end

        function obj = predict_step(obj,u,Ts)
            %{
                Performs all the steps for a priori estimation.
            %}
            
            % state prediction

            x_new = obj.predict_state(u,Ts);

            % covariance prediction

            [~,Qs,w] = obj.set_additive_noise(Ts);

            F = obj.predict_jacobian_preintegrated(u, Ts);
            G = obj.predict_covariance_preintegrated(w);

            P_new = F*obj.P*(F')+G*Qs*(G');
            P_new = 0.5*(P_new+P_new'); % WTF is this supposed to do ?

            obj.x = x_new;
            obj.P = P_new;

        end

        function obj = update_step(obj,z,h_x,H,R)
            %{
                Turns the a priori prediction into the a posteriori one for a given measurement point
            %}
            
            nx = size(obj.x,1);
            inov = z-h_x;
            S = H*obj.P*(H')+R;
            K = obj.P*(H')*inv(S);
            x_new = obj.x + K*inov;
            P_new = (eye(nx)-K*H)*obj.P;
            % NB : Not using (I-K*H)*obj.P*(I-K*H)' + K*R*K'
            % as it is implicitely included when computing S
        end
   
        function obj = update_step_sensors(obj, z)
            %{
                "Meta-function", use this to call update
                Kept separate to split math and sensors.
            %}
            R = 1 * obj.baro_measurement_uncertainty;
            h_x = obj.measurement_model();
            H = obj.measurement_jacobian();
            
            obj = obj.update_step(z, h_x, H, R);
        end

        function h_x = measurement_model(obj)
            %{
                Measurement model, currently : height (pos down)
            %}
            pd = obj.x(3);
            h_x = pd;
        end

        function H = measurement_jacobian(obj)
            %{
                Jacobian matrix of the measurement model
            %}
            H = [...
                0, 0, 1,    0, 0, 0,    0, 0, 0,    0];
        end

        function obj = EKF_rocket(x_init,init_process_cov)
            obj.x = x_init;
            obj.P = ones(10)*init_process_cov;
        end

        function [Qs, w] = generate_noise(obj)
            Fs = 1/obj.dt;

            obj.scale_var = 0.5*(1./(Fs.^2));
            obj.vel_delta_bias_sigma = obj.scale_var.* obj.accel_bias_noise;
            obj.pos_delta_bias_sigma = obj.scale_var.* obj.baro_bias_noise;

            Qs = diag([obj.additiveNoise.*ones(1,6), obj.vel_delta_bias_sigma*ones(1,3), obj.pos_delta_bias_sigma*ones(1,1)]);
            w = obj.scale_var.*[obj.accel_noise*ones(1,3), obj.baro_noise*ones(1,1)];
        end


    end

end