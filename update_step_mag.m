function   [x_new,P_new] = update_step_mag(x,P,z)

    Rmag = 0.09;
    R_mag = eye(3)*Rmag;
    
    h_x = mag_measurement_model(x);

    H_mag = mag_measurement_jacobian(x);
    
    [x_new,P_new]  = update_step(x,P,z,h_x,H_mag,R_mag);
end

function H = mag_measurement_jacobian(obj)
    q0 = obj.x(1);
    q1 = obj.x(2);
    q2 = obj.x(3);
    q3 = obj.x(4);
    magNavX = obj.x(17);
    magNavY = obj.x(18);
    magNavZ = obj.x(19);

    H = [ ...
                2*magNavY*q3 - 2*magNavZ*q2 + 2*magNavX*q0, 2*magNavZ*q3 + 2*magNavY*q2 + 2*magNavX*q1, 2*magNavY*q1 - 2*magNavZ*q0 - 2*magNavX*q2, 2*magNavZ*q1 + 2*magNavY*q0 - 2*magNavX*q3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, q0^2 + q1^2 - q2^2 - q3^2,         2*q0*q3 + 2*q1*q2,         2*q1*q3 - 2*q0*q2, 1, 0, 0;
                2*magNavZ*q1 + 2*magNavY*q0 - 2*magNavX*q3, 2*magNavZ*q0 - 2*magNavY*q1 + 2*magNavX*q2, 2*magNavZ*q3 + 2*magNavY*q2 + 2*magNavX*q1, 2*magNavZ*q2 - 2*magNavY*q3 - 2*magNavX*q0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         2*q1*q2 - 2*q0*q3, q0^2 - q1^2 + q2^2 - q3^2,         2*q0*q1 + 2*q2*q3, 0, 1, 0;
                2*magNavZ*q0 - 2*magNavY*q1 + 2*magNavX*q2, 2*magNavX*q3 - 2*magNavY*q0 - 2*magNavZ*q1, 2*magNavY*q3 - 2*magNavZ*q2 + 2*magNavX*q0, 2*magNavZ*q3 + 2*magNavY*q2 + 2*magNavX*q1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         2*q0*q2 + 2*q1*q3,         2*q2*q3 - 2*q0*q1, q0^2 - q1^2 - q2^2 + q3^2, 0, 0, 1];


   end


        function z =  mag_measurement_model(obj)
            q0 = obj.x(1);
            q1 = obj.x(2);
            q2 = obj.x(3);
            q3 = obj.x(4);
            magNavX = obj.x(17);
            magNavY = obj.x(18);
            magNavZ = obj.x(19);
            magBiasX = obj.x(20);
            magBiasY = obj.x(21);
            magBiasZ = obj.x(22);

            mx = magBiasX + magNavX*(q0^2 + q1^2 - q2^2 - q3^2) - magNavZ*(2*q0*q2 - 2*q1*q3) + magNavY*(2*q0*q3 + 2*q1*q2);
            my = magBiasY + magNavY*(q0^2 - q1^2 + q2^2 - q3^2) + magNavZ*(2*q0*q1 + 2*q2*q3) - magNavX*(2*q0*q3 - 2*q1*q2);
            mz = magBiasZ + magNavZ*(q0^2 - q1^2 - q2^2 + q3^2) - magNavY*(2*q0*q1 - 2*q2*q3) + magNavX*(2*q0*q2 + 2*q1*q3);

            z = [mx my mz]';
        end