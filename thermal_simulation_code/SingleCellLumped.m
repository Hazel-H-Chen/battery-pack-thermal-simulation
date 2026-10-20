clear; clc;

%% Load LapSim data
data   = readmatrix('lapsim_results_end_full_braking.csv');
t_data = data(:, 1);            % time (s)
P_data = data(:, 3);            % pack power (kW)
P_data(P_data < 0) = 0;         % ignores regen

%% Cell + pack parameters
m          =  0.069;      % cell mass (kg)
cp         =  900;        % specific heat (J/kg-K)
cell_data  =  readmatrix('P50B_DCIR.csv');
soc_table  =  cell_data(:, 1);
dcir_table =  cell_data(:, 2)/1000;
v_curve   = readmatrix("P50B_voltage_vs_SoC.csv");
V_nom      =  3.7;       % nominal voltage
V_max      =  4.2;       % max voltage
C_mAh      =  4600;     % cell capacity (mAh)
Q_cell     =  (C_mAh / 1000) * 3600; % calculates cell capacity in Ah

n_series  = 140;
n_parallel = 3;
n_cells   = n_parallel *n_series;        % cells in pack
hA        = 0.2;

%% Ambient / initial conditions
T_amb = 25;             % ambient temperature (C)
T0    = 30;             % initial cell temperature (C)

%% Build a multi-lap current profile (20 laps back to back)
lap_duration = t_data(end);
t_lap = t_data;
P_lap = P_data;
t_data = [];
P_data = [];
for i = 0:19
    t_data = [t_data; t_lap + i * lap_duration];
    P_data = [P_data; P_lap];
end

tspan = [t_data(1), t_data(end)];   % full 20-lap run

state0 = [T0; 1.0];

%% Solve
opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);
[t, state] = ode45(@(t, s) cell_thermal(t, s, m, cp, ...   
                             hA, T_amb, t_data, P_data, n_cells, Q_cell, soc_table, dcir_table,v_curve), tspan, state0, opts);

T   = state(:,1);
SOC = state(:,2);

%% Plot per cell Power
figure;
P_pack = interp1(t_data, P_data, t, 'linear', 'extrap') * 1000; % W
P_cell = P_pack / n_cells;
plot(t/60, P_cell, 'b', 'LineWidth',2);
xlabel('Time (min)','FontSize',20);
ylabel('Power (W)','FontSize',20);
xlim([0 0.1]);
title('Per Cell Power Per Lap','FontSize',20);
grid on;

%% Plot
figure;
C_T = plot(t/60, T, 'r', 'LineWidth', 2);
xlabel('Time (min)','FontSize',20);
ylabel('Cell Temperature (\circC)','FontSize',20);
title('Single-Cell Thermal Response Over Endurance Run','FontSize',20);
grid on;
dt1 = datatip(C_T,5,interp1(t,T,5*60),Location="southeast");
dt2 = datatip(C_T,10,interp1(t,T,10*60),Location="southeast");
dt3 = datatip(C_T,15,interp1(t,T,15*60),Location="southeast");
dt4 = datatip(C_T,20,interp1(t,T,20*60),Location="southeast");
finalTemp = datatip(C_T,T(end),interp1(t,T,20*60),Location = "southeast");

C_T.DataTipTemplate.DataTipRows(1).Label = 'Time (min):';
C_T.DataTipTemplate.DataTipRows(2).Label = 'Temperature (°C):';
C_T.DataTipTemplate.DataTipRows(1).Format = '%.1f';
C_T.DataTipTemplate.DataTipRows(2).Format = '%.1f';

fprintf('Peak cell temperature: %.1f deg C\n', max(T));
fprintf('Final cell temperature: %.1f deg C\n', T(end));

%% Plot SOC over run
figure;
SOC_plot = plot(t/60, SOC, 'g', 'LineWidth', 2);
xlabel('Time (min)','FontSize',20);
ylabel('State of Charge (%)','FontSize',20);
title('State of Charge Over Endurance Run','FontSize',20);
grid on;

%ODE Function
function dstate = cell_thermal(t, s, m, cp, ...    %cell thermal and SOC estimation
                             hA, T_amb, t_data, P_data, n_cells, Q_cell, soc_table, dcir_table,v_curve)
    T = s(1);
    soc = s(2);
    R_int = interp1(soc_table,dcir_table,soc*100,"linear","extrap");
    % Per-cell power from pack power, then current from P = V*I
    P_pack = interp1(t_data, P_data, t, 'linear', 'extrap') * 1000; % W
    P_cell = P_pack / n_cells;
    v_oc = calcV(soc,v_curve);
    v_act = calcVdrop(v_oc,P_cell,R_int);
    I      = abs(P_cell) / v_act;

    % Heat balance
    Q_gen  = I^2 * R_int;          % ohmic heat generation (W)
    Q_loss = hA * (T - T_amb);     % convective loss to ambient (W)

    dTdt = (Q_gen - Q_loss) / (m * cp);
    dSOCdt = -I/Q_cell;
    dstate = [dTdt; dSOCdt];
end

%Voltage Drop Function
function v_actual = calcVdrop(V_oc,P,IR)
v_actual = (sqrt((V_oc^2)-4*P*IR)+V_oc)/2;
end

%Get_Voltage
function v_oc = calcV(soc,v_curve)
soc_curve = v_curve(:, 1);
voltage_curve = v_curve(:, 2);
v_oc = interp1(soc_curve, voltage_curve, soc * 100, 'linear', 'extrap');
end