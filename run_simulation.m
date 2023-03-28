clc,close all,clear all
% EKF code - Samuel Wahba
% state variable: x = [ posNED, velNED, accBias ] (9x1)

%use realtime plotter from "Pose Estimation From Asynchronous Sensors" from Matlab's exemple
realtime_plots = false;

% custom plots
plot_pos_vel = true;
plot_bias = true;

%% Load data
%ld = load("drone_traj.mat");

[T, press_array, accel_array, ~] = import_data("static_test_1.csv");
press_array = table2array(press_array);
accel_array = table2array(accel_array);
data_points = size(T, 1);
%% EKF parameters
% sampling rate
Fs = 160; % FIXME Sensor
Ts = 1/Fs;

%% initial state

Nav = 100;
initstate = zeros(10,1);
P0 = ones(10)*1e-9;

ekf = EKF_rocket(initstate,P0);

%% Plotting stuff (Thanks matlab)

if(realtime_plots)
    useErrScope = true; % Turn on the streaming error plot.
    usePoseView = true; % Turn on the 3D pose viewer.
    if usePoseView
        posescope = PoseViewerWithSwitches(...
            'XPositionLimits', [-30 30], ...
            'YPositionLimits', [-30, 30], ...
            'ZPositionLimits', [-10 10]);
    end
    f = gcf;

    if useErrScope
        errscope = HelperScrollingPlotter(...
            'NumInputs', 3, ...
            'TimeSpan', 10, ...
            'SampleRate', Fs, ...
            'YLabel', {'meters', ...
            'meters', ...
            'meters'}, ...
            'Title', {'Position X Error', ...
            'Position Y Error', ...
            'Position Z Error'}, ...
            'YLimits', ...
            [-2, 2
            -2 2
            -2 2]);
    end

end


%% Simulation loop

x_traj = zeros(data_points, 10);

if(~realtime_plots)
    tic
end

for k = 1:data_points
    
    % extracting next datapoint
    press = press_array(k, :);
    vel = [0, 0, 0];
    accel = accel_array(k,:); % FIXME: '-' was there before
    
    % updating state
    ekf = ekf.predict_step(accel',Ts);
    ekf = ekf.update_step_sensors([press]');

    x_traj(k,:) = (ekf.get_state());
end

if(~realtime_plots)
    toc
end

%% Plots

% time vector
tt = (0:(size(ld,1)-1))/Fs;

if(plot_pos_vel)

    figure
    % pos
    subplot(6,1,1)
    % plot(tt,gps_traj(:,1),'--')
    plot(tt,x_traj(:,5))
    hold on
    plot(tt,ld.trajData.Position(:,1),'--')
    title("Position x axis")
    
    subplot(6,1,2)
    plot(tt,x_traj(:,6))
    hold on
    plot(tt,ld.trajData.Position(:,2),'--')
    title("Position y axis")
    
    subplot(6,1,3)
    hold on
    plot(tt,x_traj(:,7))
    plot(tt,ld.trajData.Position(:,3),'--')
    title("Position z axis")
    
    % Velocity 
    subplot(6,1,4)
    
    % plot(tt,gps_traj(:,4),'--')
    hold on
    plot(tt,x_traj(:,8))
    plot(tt,ld.trajData.Velocity(:,1),'--')
    % legend("GPS","EKF","Truth")
    title("Velocity inertial x axis")
    
    subplot(6,1,5)
    % plot(tt,gps_traj(:,5),'--')
    hold on
    plot(tt,x_traj(:,9))
    plot(tt,ld.trajData.Velocity(:,2),'--')
    % legend("GPS","EKF","Truth")
    title("Velocity inertial y axis")
    
    subplot(6,1,6)
    % plot(tt,gps_traj(:,6),'--')
    hold on
    plot(tt,x_traj(:,10))
    plot(tt,ld.trajData.Velocity(:,3),'--')
    % legend("GPS","EKF","Truth")
    title("Velocity inertial z axis")

end

if(plot_bias)
    
    % bias plots
    figure
    % magnetometer bias
    subplot(9,1,1)
    plot(tt,x_traj(:,20))
    title("Magnetometer bias x axis")
    
    subplot(9,1,2)
    plot(tt,x_traj(:,21))
    title("Magnetometer bias y axis")
    
    subplot(9,1,3)
    plot(tt,x_traj(:,22))
    title("Magnetometer bias z axis")
    
    % Gyroscope bias 
    subplot(9,1,4)
    plot(tt,x_traj(:,11))
    title("angle bias x axis")
    
    subplot(9,1,5)
    plot(tt,x_traj(:,12))
    title("angle bias y axis")
    
    subplot(9,1,6)
    plot(tt,x_traj(:,13))
    title("angle bias z axis")
    
    % Accelerometer bias
    subplot(9,1,7)
    plot(tt,x_traj(:,14))
    title("Accelerometer bias x axis")
    
    subplot(9,1,8)
    plot(tt,x_traj(:,15))
    title("Accelerometer bias y axis")
    
    subplot(9,1,9)
    plot(tt,x_traj(:,16))
    title("Accelerometer inertial z axis")

end

