function [pos, vel, accel, gyro, quat] = import_data(file)    
    T = readtable(file);
    pos = T(:, [1, 2, 3]);
    vel = T(:, [4, 5, 6]);
    accel = T(:, [47, 8, 9]);
    gyro = T(:, [10, 11, 12]);
    quat = quaternion(T(:, [13, 14, 15, 16]));
end