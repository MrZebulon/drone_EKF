function trajData = import_data(file)
    headers = ["accel_x", "accel_y", "accel_z", ...
            "gyro_x", "gyro_y", "gyro_z", ...
            "orientation"];
    
    T = readtable(file);
    as_quat = quaternion(T(:, [7, 8, 9, 10]));
    T = table([T{:, 1:6} as_quat], "VariableNames", headers);
    trajData = table2struct(T, "ToScalar",true);
end