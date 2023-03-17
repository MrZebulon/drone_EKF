function trajData = import_data(file)
    headers = ["accel", "gyro", "orientation"];
    
    T = readtable(file);
    accel = T(:, [1, 2, 3]);
    gyro = T(:, [4, 5, 6]);
    as_quat = quaternion(T(:, [7, 8, 9, 10]));
    T = table([T{:, 1:6} as_quat], "VariableNames", headers);
    trajData = table2struct(T, "ToScalar", true);
end