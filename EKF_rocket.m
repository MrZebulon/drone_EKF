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
        GyroscopeNoise = 0;

        % extra additive noise
        AccelerometerBiasNoise =  2e-4;
        GyroscopeBiasNoise = 0
        additiveNoise = 1e-8;

        scale_var = -1;
        ang_delta_bias_sigma = -1;
        vel_delta_bias_sigma = -1;

    end

    methods

        function obj = set_state(obj,x)
            obj.x = x;
        end

        function rot_mat = get_rotmat_body_2_inertial(obj)
            %{
                Returns a rotation matrix body-frame => intertial-frame
            %}

            q0 = obj.x(1);
            q1 = obj.x(2);
            q2 = obj.x(3);
            q3 = obj.x(4);

            rot_mat = zeros(3);
            rot_mat(1,1) = q0*q0+q1*q1-q2*q2-q3*q3;
            rot_mat(1,2) = 2*(q1*q2-q0*q3);
            rot_mat(1,3) = 2*(q1*q3+q0*q2);
        
            rot_mat(2,1) = 2*(q1*q2+q0*q3);
            rot_mat(2,2) = q0*q0-q1*q1+q2*q2-q3*q3;
            rot_mat(2,3) = 2*(q2*q3-q0*q1);
            
            rot_mat(3,1) = 2*(q1*q3-q0*q2);
            rot_mat(3,2) = 2*(q2*q3+q0*q1);
            rot_mat(3,3) = q0*q0-q1*q1-q2*q2+q3*q3;
        end

        function obj = set_process_cov(obj,P)
            obj.P = P;
        end

        function x = get_state(obj)
            x = obj.x;
        end

        function q = get_quat(obj)
            q = obj.x(1:4);
        end

        function p = get_posNED(obj)
            p = obj.x(5:7);
        end

        function v = get_velNED(obj)
            v = obj.x(8:10);
        end

        function omega_b = get_angular_velocity_body_frame(obj)
            v = obj.x(11:13);
        end

        function euler = get_eulerZYX(obj)
            q = obj.x(1:4);
            euler = quat2eul(q','ZYX');
        end

        function P = get_process_cov(obj)
            P = obj.P;
        end

        function F = predict_jacobian(obj,ang_delta,vel_delta,dt)
            %{
            Returns the F matrix (n by n, n = len(x)) for EKF calculations.
            
            F = jacobian w/ respect to x of function f s.t.
            x^{.} = f(x^{^}, u) + w
            
            Note that the content of this matrix has been pre-determined,
            Only numerical values are set at run-time.
            %}

            F = [...
                ];


        end

        function G = predict_process_noise(obj,w)
            %{
            Returns the G matrix (n by n, n = len(x) for EKF calculations.
            
            G = jacobian w/ respect to w of function f s.t.
            x^{.} = f(x^{^}, u) + w
            
            Note that the content of this matrix has been pre-determined,
            Only numerical values are set at run-time.
            %}

            G = [...
                ];
        end

        function x = predict_state(obj,ang_delta,vel_delta,dt)
            %{
                Realises prediction computations.
                Returns the a priori prediction. 
            %}

            x = [
                ];

            %{
                KEEP IN CASE WE MEASURE ORIENTATION

                qinit = quaternion(q0,q1,q2,q3);
                % x(1:4) = compact(normalize(qinit * quaternion(ang_delta - [dax_b, day_b, daz_b], 'rotvec')));
                delta_q = [1;(ang_delta'-[dax_b; day_b; daz_b])/2];
                x(1:4) = obj.mult_quat([q0;q1;q2;q3],delta_q);
            %}
        end

        function qn = mult_quat(obj,q1,q2)
            qnew_0 = q1(1)*q2(1)     -q1(2)*q2(2) - q1(3)*q2(3)-q1(4)*q2(4);
            qnew_1 = q1(1)*q2(2)      +q2(1)*q1(2)+ q1(3)*q2(4)-q2(3)*q1(4);
            qnew_2 = q1(1)*q2(3)      +q2(1)*q1(3) - q1(2)*q2(4)+q2(2)*q1(4);
            qnew_3 = q1(1)*q2(4)+      q2(1)*q1(4) + q1(2)*q2(3)-q2(2)*q1(3);
            qn = [qnew_0;qnew_1;qnew_2;qnew_3];
        end

        function obj = predict_step(obj,accB,omegaB,Ts)
            %{
                Realises all the steps within one cycle.
            %}
            
            % state prediction

            ang_delta= omegaB'*Ts; % reset sur l'erreur : OK
            vel_delta = accB'*Ts;

            x_new = obj.predict_state(ang_delta,vel_delta,Ts);

            % covariance prediction

            [~,Qs,w] = obj.set_additive_noise(Ts);

            G = obj.predict_process_noise(w);
            F = obj.predict_jacobian(ang_delta,vel_delta,Ts);

            P_new = F*obj.P*(F')+G+Qs;

            % FIXME : adapt here if q = 0 at all times
            x_new = obj.quaternion_normalisation(x_new);
            P_new = 0.5*(P_new+P_new');

            obj.x = x_new;
            obj.P = P_new;

        end

        function x_next = quaternion_normalisation(obj,x)
            x_next = x;
            x_next(1:4) = x(1:4)/norm(x(1:4));
        end

        function [x_new,P_new]  = update_step(obj,z,h_x,H,R)
            %{
                Truns the a priori prediction into the a posteriori one for a given measurement point
            %}
            
            nx = size(obj.x,1);
            inov = z-h_x;
            S = H*obj.P*(H')+R;
            K = obj.P*(H')*inv(S);
            x_new = obj.x + K*inov;
            P_new = (eye(nx)-K*H)*obj.P;
        end

        function q = eulZYX2quat(obj,euler)
            q = (eul2quat(euler','ZYX'))';
        end
       
        function obj = EKF_rocket(x_init,init_process_cov)
            obj.x = x_init;
            obj.P = ones(22)*init_process_cov;
        end

        function [obj,Qs,w] = set_additive_noise(obj,Ts)
            Fs = 1/Ts;

            obj.scale_var = 0.5*(1./(Fs.^2));
            obj.ang_delta_bias_sigma = obj.scale_var .* obj.GyroscopeBiasNoise;
            obj.vel_delta_bias_sigma = obj.scale_var .* obj.AccelerometerBiasNoise;

            w = obj.scale_var.*[obj.GyroscopeNoise*ones(1,3), obj.AccelerometerNoise*ones(1,3)];

            Qs = diag([obj.additiveNoise.*ones(1,10), obj.ang_delta_bias_sigma*ones(1,3), obj.vel_delta_bias_sigma*ones(1,3),  obj.GeomagneticVectorNoise*ones(1,3), obj.MagnetometerBiasNoise*ones(1,3)]);

        end


    end

end