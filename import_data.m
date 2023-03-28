function [T, press, accel, gyro] = import_data(file)    
    T = readtable(file);
    press = T(:, 2);
    accel = T(:, [3, 4, 5]);
    gyro = T(:, [6, 7, 8]);
end