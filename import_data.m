function trajData = import_data(file)
    headers = ["vel", "accel", "gyro", "orientation"];
    
    T = readtable(file);
    vel = T(:, [1, 2, 3]);
    accel = T(:, [4, 5, 6]);
    gyro = T(:, [7, 8, 9]);
    quat = quaternion(T(:, [10, 11, 12, 13]));
    T = table([vel accel gyro quat], "VariableNames", headers);
    trajData = table2struct(T, "ToScalar", true);
end