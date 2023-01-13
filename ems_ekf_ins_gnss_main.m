%   main module

%   This is a basic INS/GNSS loosely-coupled system using EKF
%   This is only a DEMO with basic updates (position/velocity) are applied
%   More advanced updates such as nonholonomic constraints, zero-speed, and
%   adaptive EKF are NOT implemented in this DEMO. The purpose of the DEMO 
%   is to demonstrate the basics of EKF in a basic INS/GNSS fusion scenario
%   For commercial use or embdded C/C++ versions, please contact mohamed.atia@carleton.ca

%   Copyright (C) 2018, Mohamed Atia, all rights reserved.
%   The software is given under GNU Lesser General Public License
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this program. If not, see
%   <http://www.gnu.org/licenses/>.

%   For commercial use or embdded C/C++ versions, please contact mohamed.atia@carleton.ca

set(0,'DefaultFigureWindowStyle','docked');
opengl('save', 'software')
clear;
close all;
fclose all;
clc;

%--> Call symbolic mechanization equation definition script
ems_symbolic_engine;

%--> Enable/Disable animation of results
enable_animation = 1;

%--> Choose the data you are ruuning (simulated or real data)
% default value will be assigned to "sim" (simulated data), change to "real" if
% real data is used.
data_type = 'sim';

%--> Converting symbolic functions into real-valued functions
C_LB_from_Euler_fun         = matlabFunction(C_LB_from_Euler);
w_L_IL_fun                  = matlabFunction(w_L_IL);
a_dot_fun                   = matlabFunction(a_dot);
b_dot_fun                   = matlabFunction(b_dot);
c_dot_fun                   = matlabFunction(c_dot);
d_dot_fun                   = matlabFunction(d_dot);
C_LB_from_quat_fun          = matlabFunction(C_LB_from_quat);
C_EN_fun                    = matlabFunction(C_EN);
C_EL_fun                    = matlabFunction(C_EL);
w_N_EN_fun                  = matlabFunction(w_N_EN);
w_N_IE_fun                  = matlabFunction(w_N_IE);
V_N_dot_fun                 = matlabFunction(V_N_dot);
pos_quat_dot_fun_EN         = matlabFunction(pos_quat_dot_EN);
pos_quat_dot_fun_EL         = matlabFunction(pos_quat_dot_EL);
C_EN_from_pos_quat_fun      = matlabFunction(C_EN_from_pos_quat);
C_EL_from_pos_quat_fun      = matlabFunction(C_EL_from_pos_quat);
pos_quat_from_lat_lon_fun   = matlabFunction(pos_quat_from_lat_lon);
w_L_IE_fun                  = matlabFunction(w_L_IE);
gravity_fun                 = matlabFunction(g);
F_fun                       = matlabFunction(F);
G_fun                       = matlabFunction(G);

%--> Earth ellipsoid shape parameters
earth_a           = 6378137;
earth_f           = 1/298.257223563;
earth_b           = earth_a*(1-earth_f);
earth_e2          = 1-(earth_b^2)/(earth_a^2);
we_value          = 2*pi/(24*60*60);

%--> Load and display data
if strcmp(data_type,'sim')
    load('ems_data_simulated.mat');
else
    load('ems_data_for_lc.mat');
end

%imupred = xlsread('imu_results.xlsx');
figure;
plot(ems_data.time,ems_data.gyro(:,1)*R2D,'r');hold on;plot(ems_data.time,ems_data.gyro(:,2)*R2D,'g');
plot(ems_data.time,ems_data.gyro(:,3)*R2D,'b');legend('gyro x','gyro y','gyro z');title('Raw IMU - Gyro(deg/s)');

figure;
plot(ems_data.time,ems_data.accel(:,1),'r');hold on;plot(ems_data.time,ems_data.accel(:,2),'g');
plot(ems_data.time,ems_data.accel(:,3),'b');legend('gyro x','gyro y','gyro z');title('Raw IMU - Accel(m/s)');

figure;
plot(ems_data.time,ems_data.roll*R2D,'r');hold on;plot(ems_data.time,ems_data.pitch*R2D,'g');
plot(ems_data.time,ems_data.heading*R2D,'b');legend('roll','pitch','heading');
xlabel('time(s)');ylabel('Orientation(deg)');title('Reference 3D Orientation');

figure;
plot(ems_data.time,ems_data.vel(:,1),'r');hold on;plot(ems_data.time,ems_data.vel(:,2),'g');
plot(ems_data.time,ems_data.vel(:,3),'b');legend('Vn','Ve','Vd');
xlabel('time(s)');ylabel('Velocity(m/s)');title('Reference 3D Velocity');

figure;
plot(ems_data.time,ems_data.east,'r');hold on;plot(ems_data.time,ems_data.north,'g');
plot(ems_data.time,ems_data.h,'b');legend('East(m)','North(m)','Altitude(m)');
xlabel('time(s)');ylabel('Position(m)');title('Reference 3D Position');

%--> Initialization
num_of_samples = length(ems_data.time);
sampling_period_sec = mean(diff(ems_data.time));
sampling_freq = round(1/sampling_period_sec);

%--> Initialize vectors sizes
time_vector_sec                         = zeros(num_of_samples,1);
raw_acc_x                               = zeros(num_of_samples,1);
raw_acc_y                               = zeros(num_of_samples,1);
raw_acc_z                               = zeros(num_of_samples,1);
raw_gyro_x                              = zeros(num_of_samples,1);
raw_gyro_y                              = zeros(num_of_samples,1);
raw_gyro_z                              = zeros(num_of_samples,1);
a_value                                 = zeros(num_of_samples,1);
b_value                                 = zeros(num_of_samples,1);
c_value                                 = zeros(num_of_samples,1);
d_value                                 = zeros(num_of_samples,1);
Euler_pitch_value                       = zeros(num_of_samples,1);
Euler_roll_value                        = zeros(num_of_samples,1);
Euler_heading_value                     = zeros(num_of_samples,1);
ve_value                                = zeros(num_of_samples,1);
vn_value                                = zeros(num_of_samples,1);
vu_value                                = zeros(num_of_samples,1);
a_pos_value                             = zeros(num_of_samples,1);
b_pos_value                             = zeros(num_of_samples,1);
c_pos_value                             = zeros(num_of_samples,1);
d_pos_value                             = zeros(num_of_samples,1);
lat_value                               = zeros(num_of_samples,1);
lon_value                               = zeros(num_of_samples,1);
alt_value                               = zeros(num_of_samples,1);
Rn_value                                = zeros(num_of_samples,1);
Rm_value                                = zeros(num_of_samples,1);
EP_value                                = zeros(num_of_samples,1);
NP_value                                = zeros(num_of_samples,1);

gyro_bias_x_value                       = zeros(num_of_samples,1);
gyro_bias_y_value                       = zeros(num_of_samples,1);
gyro_bias_z_value                       = zeros(num_of_samples,1);
acc_bias_x_value                        = zeros(num_of_samples,1);
acc_bias_y_value                        = zeros(num_of_samples,1);
acc_bias_z_value                        = zeros(num_of_samples,1);

gyro_sf_x_value                         = zeros(num_of_samples,1);
gyro_sf_y_value                         = zeros(num_of_samples,1);
gyro_sf_z_value                         = zeros(num_of_samples,1);
acc_sf_x_value                          = zeros(num_of_samples,1);
acc_sf_y_value                          = zeros(num_of_samples,1);
acc_sf_z_value                          = zeros(num_of_samples,1);

g_value                                 = zeros(num_of_samples,1);

gyro_rw_stdv_x_value                    = zeros(num_of_samples,1);
gyro_rw_stdv_y_value                    = zeros(num_of_samples,1);
gyro_rw_stdv_z_value                    = zeros(num_of_samples,1);
gyro_bias_gauss_markov_stdv_x_value     = zeros(num_of_samples,1);
gyro_bias_gauss_markov_stdv_y_value     = zeros(num_of_samples,1);
gyro_bias_gauss_markov_stdv_z_value     = zeros(num_of_samples,1);
gyro_bias_x_time_cnst_value             = zeros(num_of_samples,1);
gyro_bias_y_time_cnst_value             = zeros(num_of_samples,1);
gyro_bias_z_time_cnst_value             = zeros(num_of_samples,1);
gyro_sf_gauss_markov_stdv_x_value       = zeros(num_of_samples,1);
gyro_sf_gauss_markov_stdv_y_value       = zeros(num_of_samples,1);
gyro_sf_gauss_markov_stdv_z_value       = zeros(num_of_samples,1);
gyro_sf_x_time_cnst_value               = zeros(num_of_samples,1);
gyro_sf_y_time_cnst_value               = zeros(num_of_samples,1);
gyro_sf_z_time_cnst_value               = zeros(num_of_samples,1);

acc_rw_stdv_x_value                     = zeros(num_of_samples,1);
acc_rw_stdv_y_value                     = zeros(num_of_samples,1);
acc_rw_stdv_z_value                     = zeros(num_of_samples,1);
acc_bias_gauss_markov_stdv_x_value      = zeros(num_of_samples,1);
acc_bias_gauss_markov_stdv_y_value      = zeros(num_of_samples,1);
acc_bias_gauss_markov_stdv_z_value      = zeros(num_of_samples,1);
acc_bias_x_time_cnst_value              = zeros(num_of_samples,1);
acc_bias_y_time_cnst_value              = zeros(num_of_samples,1);
acc_bias_z_time_cnst_value              = zeros(num_of_samples,1);
acc_sf_gauss_markov_stdv_x_value        = zeros(num_of_samples,1);
acc_sf_gauss_markov_stdv_y_value        = zeros(num_of_samples,1);
acc_sf_gauss_markov_stdv_z_value        = zeros(num_of_samples,1);
acc_sf_x_time_cnst_value                = zeros(num_of_samples,1);
acc_sf_y_time_cnst_value                = zeros(num_of_samples,1);
acc_sf_z_time_cnst_value                = zeros(num_of_samples,1);

a_pos_rw_stdv_value                     = zeros(num_of_samples,1);
b_pos_rw_stdv_value                     = zeros(num_of_samples,1);
c_pos_rw_stdv_value                     = zeros(num_of_samples,1);
d_pos_rw_stdv_value                     = zeros(num_of_samples,1);
alt_rw_stdv_value                       = zeros(num_of_samples,1);

wg_noise_value                          = zeros(num_of_samples,1);

%--> Initialize the state
Euler_roll_value(1)     = ems_data.roll(1);
Euler_pitch_value(1)    = ems_data.pitch(1);
Euler_heading_value(1)  = ems_data.heading(1);
C_LB_value = C_LB_from_Euler_fun(Euler_roll_value(1),Euler_pitch_value(1),Euler_heading_value(1));
attitude_quat_value = convert_dcm_to_quat (C_LB_value);
a_value(1) = attitude_quat_value(1);
b_value(1) = attitude_quat_value(2);
c_value(1) = attitude_quat_value(3);
d_value(1) = attitude_quat_value(4);
lat_value(1)  = ems_data.lat(1)*R2D;
lon_value(1)  = ems_data.lon(1)*R2D;
alt_value(1)  = ems_data.h(1);
Rn_value(1)   = earth_a/sqrt(1-earth_e2*sin(lat_value(1)*D2R)*sin(lat_value(1)*D2R));
Rm_value(1)   = earth_a*(1-earth_e2)/((1-earth_e2*sin(lat_value(1)*D2R)*sin(lat_value(1)*D2R))^(1.5));
EP_value(1) = 0;
NP_value(1) = 0;
g_value(1) = gravity_fun(Rm_value(1),Rn_value(1),alt_value(1),lat_value(1));
C_EN_value = C_EN_fun(lat_value(1)*D2R , lon_value(1)*D2R);
pos_quat_vector = convert_dcm_to_quat (C_EN_value);
a_pos_value(1) = pos_quat_vector(1)/sqrt(sum(pos_quat_vector.^2));
b_pos_value(1) = pos_quat_vector(2)/sqrt(sum(pos_quat_vector.^2));
c_pos_value(1) = pos_quat_vector(3)/sqrt(sum(pos_quat_vector.^2));
d_pos_value(1) = pos_quat_vector(4)/sqrt(sum(pos_quat_vector.^2));
ve_value(1) = ems_data.vel_N(1);
vn_value(1) = ems_data.vel_N(2);
vu_value(1) = ems_data.vel_N(3);

%--> Initial biases
gyro_bias_x_value(1) = 0.0;
gyro_bias_y_value(1) = 0.0;
gyro_bias_z_value(1) = 0.0;
acc_bias_x_value(1)  = 0.0;
acc_bias_y_value(1)  = 0.0;
acc_bias_z_value(1)  = 0.0;

%--> Initial scale factors
gyro_sf_x_value(1) = 0.0;
gyro_sf_y_value(1) = 0.0;
gyro_sf_z_value(1) = 0.0;
acc_sf_x_value(1)  = 0.0;
acc_sf_y_value(1)  = 0.0;
acc_sf_z_value(1)  = 0.0;

%--> Position system noise in Q matrix
a_pos_rw_stdv_value(1)  = 0.0;
b_pos_rw_stdv_value(1)  = 0.0;
c_pos_rw_stdv_value(1)  = 0.0;
d_pos_rw_stdv_value(1)  = 0.0;
alt_rw_stdv_value(1)    = 0.0;

%--> Accelerometer System noise and Gauss-Markov model parameters
acc_rw_stdv_x_value(1) = 0.01;
acc_rw_stdv_y_value(1) = 0.01;
acc_rw_stdv_z_value(1) = 0.01;
acc_bias_gauss_markov_stdv_x_value (1) = 0.01;
acc_bias_gauss_markov_stdv_y_value (1) = 0.01;
acc_bias_gauss_markov_stdv_z_value (1) = 0.01;
acc_bias_x_time_cnst_value(1) = 3*3600;
acc_bias_y_time_cnst_value(1) = 3*3600;
acc_bias_z_time_cnst_value(1) = 3*3600;
acc_sf_gauss_markov_stdv_x_value (1) = 1e-12;
acc_sf_gauss_markov_stdv_y_value (1) = 1e-12;
acc_sf_gauss_markov_stdv_z_value (1) = 1e-12;
acc_sf_x_time_cnst_value(1) = 3*3600;
acc_sf_y_time_cnst_value(1) = 3*3600;
acc_sf_z_time_cnst_value(1) = 3*3600;

%--> Gyroscope System noise and Gauss-Markov model parameters
gyro_rw_stdv_x_value(1) = 0.01;
gyro_rw_stdv_y_value(1) = 0.01;
gyro_rw_stdv_z_value(1) = 0.01;
gyro_bias_gauss_markov_stdv_x_value (1) = 0.001*D2R;
gyro_bias_gauss_markov_stdv_y_value (1) = 0.001*D2R;
gyro_bias_gauss_markov_stdv_z_value (1) = 0.001*D2R;
gyro_bias_x_time_cnst_value(1) = 3*3600;
gyro_bias_y_time_cnst_value(1) = 3*3600;
gyro_bias_z_time_cnst_value(1) = 3*3600;
gyro_sf_gauss_markov_stdv_x_value (1) = 1e-12;
gyro_sf_gauss_markov_stdv_y_value (1) = 1e-12;
gyro_sf_gauss_markov_stdv_z_value (1) = 1e-12;
gyro_sf_x_time_cnst_value(1) = 3*3600;
gyro_sf_y_time_cnst_value(1) = 3*3600;
gyro_sf_z_time_cnst_value(1) = 3*3600;

%--> Noise is to set to zero in INS prediction
wg_noise_value(1) = 0.0;

%--> Initial error states covariances
east_pos_error_covariance(1) = 25^2;
north_pos_error_covariance(1) = 25^2;
alt_error_covariance(1) = 20^2;
ve_error_covariance(1) = 50.05^2;
vn_error_covariance(1) = 50.05^2;
vu_error_covariance(1) = 50.05^2;
b_error_covariance(1) = 2.5^2;
c_error_covariance(1) = 2.5^2;
d_error_covariance(1) = 2.5^2;
gyro_bias_x_error_covariance(1) = .2^2;
gyro_bias_y_error_covariance(1) = .4^2;
gyro_bias_z_error_covariance(1) = .6^2;
acc_bias_x_error_covariance(1) = 0.05^2;
acc_bias_y_error_covariance(1) = 0.05^2;
acc_bias_z_error_covariance(1) = 0.05^2;
gyro_sf_x_error_covariance(1) = 0.001^2;
gyro_sf_y_error_covariance(1) = 0.001^2;
gyro_sf_z_error_covariance(1) = 0.001^2;
acc_sf_x_error_covariance(1) = 0.001^2;
acc_sf_y_error_covariance(1) = 0.001^2;
acc_sf_z_error_covariance(1) = 0.001^2;

P = diag([
    east_pos_error_covariance(1);
    north_pos_error_covariance(1);
    alt_error_covariance(1);
    ve_error_covariance(1);
    vn_error_covariance(1);
    vu_error_covariance(1);
    b_error_covariance(1);
    c_error_covariance(1);
    d_error_covariance(1);
    gyro_bias_x_error_covariance(1);
    gyro_bias_y_error_covariance(1);
    gyro_bias_z_error_covariance(1);
    acc_bias_x_error_covariance(1);
    acc_bias_y_error_covariance(1);
    acc_bias_z_error_covariance(1);
    gyro_sf_x_error_covariance(1);
    gyro_sf_y_error_covariance(1);
    gyro_sf_z_error_covariance(1);
    acc_sf_x_error_covariance(1);
    acc_sf_y_error_covariance(1);
    acc_sf_z_error_covariance(1)
    ]);

%--> Draw a 3D Aircraft for animation
X = 60;Y = 15;Z = 5;
origin = [X/2 Y/2 Z/2];
initial_vert = ...
   [X Y 0;              %(1)
    0 Y 0;              %(2)
    0 Y Z;              %(3)
    X Y Z;              %(4)
    0 0 Z;              %(5)
    X 0 Z;              %(6)
    X 0 0;              %(7)
    0 0 0;              %(8)
    1.3*X Y/2 0;        %(9)
    1.3*X Y/2 0.5*Z;    %(10)
];
for p = 1:length(initial_vert)
    initial_vert(p,:) = initial_vert(p,:) - origin;
end
CubePoints = [initial_vert(:,1),initial_vert(:,2),initial_vert(:,3)];
[ CubeXData , CubeYData , CubeZData ] = get_cube_axis_data(initial_vert);
faces = [1 2 3 4; 4 3 5 6; 6 7 8 5; 1 2 8 7; 6 7 1 4; 2 3 5 8;1 9 10 4;4 10 6 6;6 10 9 7;1 9 7 7];

%--> Loop counters
line_counter = 0;
innovation_seq_counter = 1;

%--> Main processing loop
for index = 1:length(ems_data.time)-1
    
    %--> Read raw IMU
    time_vector_sec(index) = ems_data.time(index);
    raw_acc_x(index) = ems_data.accel(index,1);
    raw_acc_y(index) = ems_data.accel(index,2);
    raw_acc_z(index) = ems_data.accel(index,3);
    raw_gyro_x(index) = ems_data.gyro(index,1);
    raw_gyro_y(index) = ems_data.gyro(index,2);
    raw_gyro_z(index) = ems_data.gyro(index,3);
    
    %--> Calculate Earth Rotation Rate, Transport Rate, and Corriollis effect
    w_N_EN_value = w_N_EN_fun(Rm_value(index),Rn_value(index),alt_value(index),lat_value(index)*D2R,ve_value(index),vn_value(index));
    w_N_IE_value = w_N_IE_fun(lat_value(index)*D2R,we_value);
    w_L_IL_value = C_LN*(w_N_EN_value + w_N_IE_value);
    w_B_IL_value = C_LB_value'*w_L_IL_value;
    corriollis_effect_vector_in_B = C_LB_value'*C_LN*cross((w_N_EN_value + 2*w_N_IE_value),[ve_value(index);vn_value(index);vu_value(index)]);
    
    %--> Perform attitude quaternion mechanization
    a_dot_value = a_dot_fun(Rm_value(index),Rn_value(index),alt_value(index),b_value(index),c_value(index),d_value(index),raw_gyro_x(index),raw_gyro_y(index),raw_gyro_z(index),gyro_sf_x_value(index),gyro_sf_y_value(index),gyro_sf_z_value(index),gyro_bias_x_value(index),gyro_bias_y_value(index),gyro_bias_z_value(index),0.0,0.0,0.0,lat_value(index)*D2R,ve_value(index),vn_value(index),we_value,wg_noise_value(index));
    b_dot_value = b_dot_fun(Rm_value(index),Rn_value(index),alt_value(index),b_value(index),c_value(index),d_value(index),raw_gyro_x(index),raw_gyro_y(index),raw_gyro_z(index),gyro_sf_x_value(index),gyro_sf_y_value(index),gyro_sf_z_value(index),gyro_bias_x_value(index),gyro_bias_y_value(index),gyro_bias_z_value(index),0.0,0.0,0.0,lat_value(index)*D2R,ve_value(index),vn_value(index),we_value,wg_noise_value(index));
    c_dot_value = c_dot_fun(Rm_value(index),Rn_value(index),alt_value(index),b_value(index),c_value(index),d_value(index),raw_gyro_x(index),raw_gyro_y(index),raw_gyro_z(index),gyro_sf_x_value(index),gyro_sf_y_value(index),gyro_sf_z_value(index),gyro_bias_x_value(index),gyro_bias_y_value(index),gyro_bias_z_value(index),0.0,0.0,0.0,lat_value(index)*D2R,ve_value(index),vn_value(index),we_value,wg_noise_value(index));
    d_dot_value = d_dot_fun(Rm_value(index),Rn_value(index),alt_value(index),b_value(index),c_value(index),d_value(index),raw_gyro_x(index),raw_gyro_y(index),raw_gyro_z(index),gyro_sf_x_value(index),gyro_sf_y_value(index),gyro_sf_z_value(index),gyro_bias_x_value(index),gyro_bias_y_value(index),gyro_bias_z_value(index),0.0,0.0,0.0,lat_value(index)*D2R,ve_value(index),vn_value(index),we_value,wg_noise_value(index));
    
    if (~isreal([a_dot_value b_dot_value c_dot_value d_dot_value]))
        disp('imaginary quaternion rate');return;
    end
    
    %--> Advance quaternion
    a_value(index+1) = a_value(index) + a_dot_value*sampling_period_sec;
    b_value(index+1) = b_value(index) + b_dot_value*sampling_period_sec;
    c_value(index+1) = c_value(index) + c_dot_value*sampling_period_sec;
    d_value(index+1) = d_value(index) + d_dot_value*sampling_period_sec;
    
    %--> Normalize the quaternion
    normalization_factor = sqrt(a_value(index+1)^2 + b_value(index+1)^2 + c_value(index+1)^2 + d_value(index+1)^2);
    a_value(index+1) = a_value(index+1)/normalization_factor;
    b_value(index+1) = b_value(index+1)/normalization_factor;
    c_value(index+1) = c_value(index+1)/normalization_factor;
    d_value(index+1) = d_value(index+1)/normalization_factor;
    
    %--> Calculate DCM from Quaternion
    C_LB_value = C_LB_from_quat_fun(...
        b_value(index+1),...
        c_value(index+1),...
        d_value(index+1));
    
    if (~isreal(C_LB_value))
        disp('imaginary DCM C_LB');return;
    end
    
    %--> Normalize the DCM
    C_LB_value(1,:)=C_LB_value(1,:)/norm(C_LB_value(1,:));
    C_LB_value(2,:)=C_LB_value(2,:)/norm(C_LB_value(2,:));
    C_LB_value(3,:)=C_LB_value(3,:)/norm(C_LB_value(3,:));
    
    %--> Calculate Euler angles from DCM
    Euler_pitch_value(index+1)     = atan2(-C_LB_value(3,1),(sqrt(C_LB_value(3,2)^2 + C_LB_value(3,3)^2)));
    Euler_roll_value(index+1)      = atan2(C_LB_value(3,2),C_LB_value(3,3));
    Euler_heading_value(index+1)   = atan2(C_LB_value(2,1),C_LB_value(1,1));
    
    if (~isreal([Euler_roll_value(index+1) Euler_pitch_value(index+1) Euler_heading_value(index+1)]))
        disp('imaginary Euler angles');return;
    end
    
    %--> Advance velocity
    V_N_dot_value = V_N_dot_fun(Rm_value(index),Rn_value(index),raw_acc_x(index),raw_acc_y(index),raw_acc_z(index),...
        acc_sf_x_value(index),acc_sf_y_value(index),acc_sf_z_value(index),...
        acc_bias_x_value(index),acc_bias_y_value(index),acc_bias_z_value(index),...
        0.0,0.0,0.0,...
        alt_value(index),b_value(index),c_value(index),d_value(index),lat_value(index)*D2R,ve_value(index),vn_value(index),vu_value(index),we_value,wg_noise_value(index));
    
    if (~isreal(V_N_dot_value))
        disp('imaginary velocity rate');return;
    end
    
    ve_value(index+1) = ve_value(index) + V_N_dot_value(1)*sampling_period_sec;
    vn_value(index+1) = vn_value(index) + V_N_dot_value(2)*sampling_period_sec;
    vu_value(index+1) = vu_value(index) + V_N_dot_value(3)*sampling_period_sec;
    
    %--> Advance position
    pos_quat_dot_value = pos_quat_dot_fun_EN(Rm_value(index),Rn_value(index),a_pos_rw_stdv_value(index), alt_value(index),...
        b_pos_value(index),b_pos_rw_stdv_value(index), c_pos_value(index),c_pos_rw_stdv_value(index), d_pos_value(index),d_pos_rw_stdv_value(index),...
        lat_value(index)*D2R,ve_value(index),vn_value(index), wg_noise_value(index));
    
    if (~isreal(pos_quat_dot_value))
        disp('imaginary position quaternion rate');return;
    end
    
    a_pos_value(index+1) = a_pos_value(index) + pos_quat_dot_value(1)*sampling_period_sec;
    b_pos_value(index+1) = b_pos_value(index) + pos_quat_dot_value(2)*sampling_period_sec;
    c_pos_value(index+1) = c_pos_value(index) + pos_quat_dot_value(3)*sampling_period_sec;
    d_pos_value(index+1) = d_pos_value(index) + pos_quat_dot_value(4)*sampling_period_sec;
    normalization_factor = sqrt(a_pos_value(index+1)^2 + b_pos_value(index+1)^2 + c_pos_value(index+1)^2 + d_pos_value(index+1)^2);
    a_pos_value(index+1) = a_pos_value(index+1)/normalization_factor;
    b_pos_value(index+1) = b_pos_value(index+1)/normalization_factor;
    c_pos_value(index+1) = c_pos_value(index+1)/normalization_factor;
    d_pos_value(index+1) = d_pos_value(index+1)/normalization_factor;
    
    %--> Calculate position DCM from position quaternion
    C_EN_value = C_EN_from_pos_quat_fun(b_pos_value(index+1),c_pos_value(index+1),d_pos_value(index+1));
    
    if (~isreal(C_EN_value))
        disp('imaginary position DCM');return;
    end
    
    %--> Calculate lat and lon from position DCM
    lon_from_C_EN_value = atan2(C_EN_value(1,3),C_EN_value(3,3));
    lat_from_C_EN_value = atan2(C_EN_value(2,3),sqrt(C_EN_value(2,1)^2+C_EN_value(2,2)^2));
    lat_value(index+1) = lat_from_C_EN_value*R2D;
    lon_value(index+1) = lon_from_C_EN_value*R2D;
    
    %--> Calculate alt
    alt_value(index+1) = alt_value(index) + vu_value(index)*sampling_period_sec;
    
    %--> Calculate new earth parameters and east, north position components
    Rn_value(index+1)  = earth_a/sqrt(1-earth_e2*sin(lat_value(index+1)*D2R)*sin(lat_value(index+1)*D2R));
    Rm_value(index+1)  = earth_a*(1-earth_e2)/((1-earth_e2*sin(lat_value(index+1)*D2R)*sin(lat_value(index+1)*D2R))^(1.5));
    EP_value(index+1) = (lon_value(index+1)-lon_value(1))*D2R*(Rn_value(index+1)+alt_value(index+1))*cos(lat_value(index+1)*D2R);
    NP_value(index+1) = (lat_value(index+1)-lat_value(1))*D2R*(Rm_value(index+1)+alt_value(index+1));
    
    %--> Advance biases
    gyro_bias_x_value(index+1) = gyro_bias_x_value(index);
    gyro_bias_y_value(index+1) = gyro_bias_y_value(index);
    gyro_bias_z_value(index+1) = gyro_bias_z_value(index);
    acc_bias_x_value(index+1) = acc_bias_x_value(index);
    acc_bias_y_value(index+1) = acc_bias_y_value(index);
    acc_bias_z_value(index+1) = acc_bias_z_value(index);
    
    %--> Advance scale factors
    gyro_sf_x_value(index+1) = gyro_sf_x_value(index);
    gyro_sf_y_value(index+1) = gyro_sf_y_value(index);
    gyro_sf_z_value(index+1) = gyro_sf_z_value(index);
    acc_sf_x_value(index+1) = acc_sf_x_value(index);
    acc_sf_y_value(index+1) = acc_sf_y_value(index);
    acc_sf_z_value(index+1) = acc_sf_z_value(index);
    
    %--> Kalman Filter model matrices calculation
    F_value = F_fun(Rm_value(index),Rn_value(index),raw_acc_x(index),raw_acc_y(index),raw_acc_z(index),...
        acc_sf_x_value(index),acc_sf_y_value(index),acc_sf_z_value(index),...
        acc_bias_x_value(index),acc_bias_y_value(index),acc_bias_z_value(index),...
        acc_rw_stdv_x_value(index),acc_rw_stdv_y_value(index),acc_rw_stdv_z_value(index),...
        acc_sf_x_time_cnst_value(index),acc_sf_y_time_cnst_value(index),acc_sf_z_time_cnst_value(index),...
        acc_bias_x_time_cnst_value(index),acc_bias_y_time_cnst_value(index),acc_bias_z_time_cnst_value(index),...
        alt_value(index),b_value(index),c_value(index),d_value(index),...
        raw_gyro_x(index),raw_gyro_y(index),raw_gyro_z(index),...
        gyro_sf_x_value(index),gyro_sf_y_value(index),gyro_sf_z_value(index),...
        gyro_bias_x_value(index),gyro_bias_y_value(index),gyro_bias_z_value(index),...
        gyro_rw_stdv_x_value(index),gyro_rw_stdv_y_value(index),gyro_rw_stdv_z_value(index),...
        gyro_sf_x_time_cnst_value(index),gyro_sf_y_time_cnst_value(index),gyro_sf_z_time_cnst_value(index),...
        gyro_bias_x_time_cnst_value(index),gyro_bias_y_time_cnst_value(index),gyro_bias_z_time_cnst_value(index),...
        lat_value(index)*R2D,ve_value(index),vn_value(index),vu_value(index),we_value,wg_noise_value(index));
    
    if (~isreal(F_value))
        disp('imaginary transition matrix');return;
    end
    
    Phi = eye(size(F_value)) + F_value*sampling_period_sec;
    
    G_value = G_fun(acc_rw_stdv_x_value(index),acc_rw_stdv_y_value(index),acc_rw_stdv_z_value(index),...
        acc_sf_x_time_cnst_value(index),acc_sf_y_time_cnst_value(index),acc_sf_z_time_cnst_value(index),...
        acc_bias_x_time_cnst_value(index),acc_bias_z_time_cnst_value(index),acc_bias_y_time_cnst_value(index),...
        acc_sf_gauss_markov_stdv_x_value(index),acc_sf_gauss_markov_stdv_y_value(index),acc_sf_gauss_markov_stdv_z_value(index),...
        acc_bias_gauss_markov_stdv_x_value(index),acc_bias_gauss_markov_stdv_y_value(index),acc_bias_gauss_markov_stdv_z_value(index),....
        alt_rw_stdv_value(index),...
        b_value(index), c_value(index),d_value(index),...
        gyro_rw_stdv_x_value(index),gyro_rw_stdv_y_value(index),gyro_rw_stdv_z_value(index),...
        gyro_sf_x_time_cnst_value(index),gyro_sf_y_time_cnst_value(index),gyro_sf_z_time_cnst_value(index),...
        gyro_bias_x_time_cnst_value(index),gyro_bias_y_time_cnst_value(index),gyro_bias_z_time_cnst_value(index),...
        gyro_sf_gauss_markov_stdv_x_value(index),gyro_sf_gauss_markov_stdv_y_value(index),gyro_sf_gauss_markov_stdv_z_value(index),...
        gyro_bias_gauss_markov_stdv_x_value(index),gyro_bias_gauss_markov_stdv_y_value(index),gyro_bias_gauss_markov_stdv_z_value(index));
    
    if (~isreal(G_value))
        disp('imaginary noise shaping matrix');return;
    end
    
    Qd = sampling_period_sec^2*(G_value*G_value');
    
    %--> EKF prediction
    P = Phi*P*Phi' + Qd;
    
    if (~isreal(P))
        disp('imaginary noise covariance matrix');return;
    end
    
    %--> apply updates from noisy obervations
    if (mod(index,sampling_freq) == 0)
        %--> Initialize error vector and matrices
        clear z H R;
        state_correction_vector = zeros(length(P),1);
        
        %--> set the innovation sequence       
        z(1) = ems_data.east_noisy(index)        - EP_value(index+1);
        z(2) = ems_data.north_noisy(index)       - NP_value(index+1);
        z(3) = ems_data.h_noisy(index)          - alt_value(index+1);
        z(4) = ems_data.vel_N_noisy(index,1)        - ve_value(index+1);
        z(5) = ems_data.vel_N_noisy(index,2)       - vn_value(index+1);
        z(6) = ems_data.vel_N_noisy(index,3)          - vu_value(index+1);
        
        H              = zeros(6,21);
        H(1,1)         = 1;
        H(2,2)         = 1;
        H(3,3)         = 1;
        H(4,4)         = 1;
        H(5,5)         = 1;
        H(6,6)         = 1;
        R              = diag([10 10 10 .5 .5 .5].^2);
        
        K = P*H'/(H*P*H'+R);
        state_correction_vector = state_correction_vector + K*z';
        P = P - K*H*P;
               
        if (~isreal(P))
            disp('imaginary updated error covariance matrix ');return;
        end
        
        correct_states;
    end
    
    %--> Normalize P
    P = (P+P')/2;
    
    if (~isreal(P))
        disp('imaginary error covariance matrix ');return;
    end
    
    gyro_rw_stdv_x_value(index+1) = gyro_rw_stdv_x_value(index);
    gyro_rw_stdv_y_value(index+1) = gyro_rw_stdv_y_value(index);
    gyro_rw_stdv_z_value(index+1) = gyro_rw_stdv_z_value(index);   
    gyro_bias_gauss_markov_stdv_x_value (index+1) = gyro_bias_gauss_markov_stdv_x_value (index);
    gyro_bias_gauss_markov_stdv_y_value (index+1) = gyro_bias_gauss_markov_stdv_y_value (index);
    gyro_bias_gauss_markov_stdv_z_value (index+1) = gyro_bias_gauss_markov_stdv_z_value (index);
    gyro_bias_x_time_cnst_value(index+1) = gyro_bias_x_time_cnst_value(index);
    gyro_bias_y_time_cnst_value(index+1) = gyro_bias_y_time_cnst_value(index);
    gyro_bias_z_time_cnst_value(index+1) = gyro_bias_z_time_cnst_value(index);
    gyro_sf_gauss_markov_stdv_x_value (index+1) = gyro_sf_gauss_markov_stdv_x_value (index);
    gyro_sf_gauss_markov_stdv_y_value (index+1) = gyro_sf_gauss_markov_stdv_y_value (index);
    gyro_sf_gauss_markov_stdv_z_value (index+1) = gyro_sf_gauss_markov_stdv_z_value (index);
    gyro_sf_x_time_cnst_value(index+1) = gyro_sf_x_time_cnst_value(index);
    gyro_sf_y_time_cnst_value(index+1) = gyro_sf_y_time_cnst_value(index);
    gyro_sf_z_time_cnst_value(index+1) = gyro_sf_z_time_cnst_value(index);
    
    acc_rw_stdv_x_value(index+1) = acc_rw_stdv_x_value(index);
    acc_rw_stdv_y_value(index+1) = acc_rw_stdv_y_value(index);
    acc_rw_stdv_z_value(index+1) = acc_rw_stdv_z_value(index);
    acc_bias_gauss_markov_stdv_x_value (index+1) = acc_bias_gauss_markov_stdv_x_value (index);
    acc_bias_gauss_markov_stdv_y_value (index+1) = acc_bias_gauss_markov_stdv_y_value (index);
    acc_bias_gauss_markov_stdv_z_value (index+1) = acc_bias_gauss_markov_stdv_z_value (index);
    acc_bias_x_time_cnst_value(index+1) = acc_bias_x_time_cnst_value(index);
    acc_bias_y_time_cnst_value(index+1) = acc_bias_y_time_cnst_value(index);
    acc_bias_z_time_cnst_value(index+1) = acc_bias_z_time_cnst_value(index);
    acc_sf_gauss_markov_stdv_x_value (index+1) = acc_sf_gauss_markov_stdv_x_value (index);
    acc_sf_gauss_markov_stdv_y_value (index+1) = acc_sf_gauss_markov_stdv_y_value (index);
    acc_sf_gauss_markov_stdv_z_value (index+1) = acc_sf_gauss_markov_stdv_z_value (index);
    acc_sf_x_time_cnst_value(index+1) = acc_sf_x_time_cnst_value(index);
    acc_sf_y_time_cnst_value(index+1) = acc_sf_y_time_cnst_value(index);
    acc_sf_z_time_cnst_value(index+1) = acc_sf_z_time_cnst_value(index);
        
    a_pos_rw_stdv_value(index + 1) = a_pos_rw_stdv_value(index);
    b_pos_rw_stdv_value(index + 1) = b_pos_rw_stdv_value(index);
    c_pos_rw_stdv_value(index + 1) = c_pos_rw_stdv_value(index);
    d_pos_rw_stdv_value(index + 1) = d_pos_rw_stdv_value(index);
    alt_rw_stdv_value(index + 1)   = alt_rw_stdv_value(index);
    
    wg_noise_value(index+1) = 0.0;

    g_value(index + 1) = gravity_fun(Rm_value(index),Rn_value(index),alt_value(index),lat_value(index));
    
    %--> Disply results in 3D
    
    C_LB_value_plus_90 = C_LB_from_Euler_fun(Euler_roll_value(index),...
        -Euler_pitch_value(index),...
        pi/2 - Euler_heading_value(index));
    
    updated_vert = C_LB_value_plus_90*initial_vert';
    
    updated_vert = updated_vert';
    
    for p = 1:length(updated_vert)
        updated_vert(p,:) = updated_vert(p,:) + [EP_value(index), NP_value(index), alt_value(index)];
    end
    
    [ CubeXData , CubeYData , CubeZData ] = get_cube_axis_data(updated_vert);
    CubePoints = [updated_vert(:,1),updated_vert(:,2),updated_vert(:,3)];
    
    %--> Create first figure after the 5th IMU sample
    if(index == 5)
        figure;
        h_heading = plot(time_vector_sec(1:index)-time_vector_sec (1), Euler_heading_value(1:index)*R2D, 'r');
        legend('heading');grid on;xlabel('time(sec)'); ylabel('heading(deg)');
        set(h_heading,'YDataSource','Euler_heading_value(1:index)*R2D');
        set(h_heading,'XDataSource','time_vector_sec(1:index)-time_vector_sec (1)');
        
        figure;
        h_roll = plot(time_vector_sec(1:index)-time_vector_sec (1), Euler_roll_value(1:index)*R2D, 'b');hold on;
        h_pitch = plot(time_vector_sec(1:index)-time_vector_sec (1), Euler_pitch_value(1:index)*R2D, 'r');
        legend('roll','pitch');grid on;xlabel('time(sec)'); ylabel('roll/pitch(deg)');
        set(h_pitch,'YDataSource','Euler_pitch_value(1:index)*R2D');
        set(h_pitch,'XDataSource','time_vector_sec(1:index)-time_vector_sec (1)');
        set(h_roll,'YDataSource','Euler_roll_value(1:index)*R2D');
        set(h_roll,'XDataSource','time_vector_sec(1:index)-time_vector_sec (1)');
               
        figure;
        h_Ve = plot(time_vector_sec(1:index)-time_vector_sec (1), ve_value(1:index), 'b');hold on;
        h_Vn = plot(time_vector_sec(1:index)-time_vector_sec (1), vn_value(1:index), 'g');
        h_Vu = plot(time_vector_sec(1:index)-time_vector_sec (1), vu_value(1:index), 'r');
        legend('ve','vn','vu');grid on;xlabel('time(sec)'); ylabel('Velocities(m/s)');
        set(h_Ve,'YDataSource','ve_value(1:index)');
        set(h_Ve,'XDataSource','time_vector_sec(1:index)-time_vector_sec (1)');
        set(h_Vn,'YDataSource','vn_value(1:index)');
        set(h_Vn,'XDataSource','time_vector_sec(1:index)-time_vector_sec (1)');
        set(h_Vu,'YDataSource','vu_value(1:index)');
        set(h_Vu,'XDataSource','time_vector_sec(1:index)-time_vector_sec (1)');
        
        figure;
        h_east = plot(time_vector_sec(1:index)-time_vector_sec (1), EP_value(1:index), 'b');hold on;
        h_north = plot(time_vector_sec(1:index)-time_vector_sec (1), NP_value(1:index), 'g');
        h_alt = plot(time_vector_sec(1:index)-time_vector_sec (1), alt_value(1:index), 'r');
        legend('east','north','up');grid on;xlabel('time(sec)'); ylabel('Position ENU(m)');
        set(h_east,'YDataSource','EP_value(1:index)');
        set(h_east,'XDataSource','time_vector_sec(1:index)-time_vector_sec (1)');
        set(h_north,'YDataSource','NP_value(1:index)');
        set(h_north,'XDataSource','time_vector_sec(1:index)-time_vector_sec (1)');
        set(h_alt,'YDataSource','alt_value(1:index)');
        set(h_alt,'XDataSource','time_vector_sec(1:index)-time_vector_sec (1)');
        
        figure;
        h_gyro_bias_x = plot(time_vector_sec(1:index)-time_vector_sec (1), gyro_bias_x_value(1:index)*R2D, 'r');grid on;hold on;
        h_gyro_bias_y = plot(time_vector_sec(1:index)-time_vector_sec (1), gyro_bias_y_value(1:index)*R2D, 'g');grid on;hold on;
        h_gyro_bias_z = plot(time_vector_sec(1:index)-time_vector_sec (1), gyro_bias_z_value(1:index)*R2D, 'b');grid on;hold on;
        legend('gyr bias x','gyr bias y', 'gyr bias z');grid on;xlabel('time(sec)'); ylabel('gyro bias (deg/s)');
        set(h_gyro_bias_x,'YDataSource','gyro_bias_x_value(1:index)*R2D');
        set(h_gyro_bias_x,'XDataSource','time_vector_sec(1:index)-time_vector_sec (1)');
        set(h_gyro_bias_y,'YDataSource','gyro_bias_y_value(1:index)*R2D');
        set(h_gyro_bias_y,'XDataSource','time_vector_sec(1:index)-time_vector_sec (1)');
        set(h_gyro_bias_z,'YDataSource','gyro_bias_z_value(1:index)*R2D');
        set(h_gyro_bias_z,'XDataSource','time_vector_sec(1:index)-time_vector_sec (1)');
        
        
        figure;
        h_accel_bias_x = plot(time_vector_sec(1:index)-time_vector_sec (1), acc_bias_x_value(1:index), 'r');grid on;hold on;
        h_accel_bias_y = plot(time_vector_sec(1:index)-time_vector_sec (1), acc_bias_y_value(1:index), 'g');grid on;hold on;
        h_accel_bias_z = plot(time_vector_sec(1:index)-time_vector_sec (1), acc_bias_z_value(1:index), 'b');grid on;hold on;
        legend('acc bias x','acc bias y', 'acc bias z');grid on;xlabel('time(sec)'); ylabel('acc bias (m/s2)');
        set(h_accel_bias_x,'YDataSource','acc_bias_x_value(1:index)');
        set(h_accel_bias_x,'XDataSource','time_vector_sec(1:index)-time_vector_sec (1)');
        set(h_accel_bias_y,'YDataSource','acc_bias_y_value(1:index)');
        set(h_accel_bias_y,'XDataSource','time_vector_sec(1:index)-time_vector_sec (1)');
        set(h_accel_bias_z,'YDataSource','acc_bias_z_value(1:index)');
        set(h_accel_bias_z,'XDataSource','time_vector_sec(1:index)-time_vector_sec (1)');
        if strcmp(data_type,'sim')
            axis([0 500 -0.5 0.8]);
        end

                        
        figure;
        hold on;
        plot3(ems_data.east,ems_data.north,ems_data.h, 'LineWidth', 1, 'color', [0.93 .69 .13]);grid on;xlabel('east');ylabel('north');zlabel('alt');
        if strcmp(data_type,'sim')
            view(-26,69);
        else
            view(20,60);
        end
        title('3D Trajectory shoing Ground Truth, Noisy Updates, and EKF output');
        h_noisy_updates = plot3(ems_data.east_noisy(1:20:index),ems_data.north_noisy(1:20:index),ems_data.h_noisy(1:20:index), '-*','MarkerSize', 5,'LineWidth', 0.5, 'color','k');grid on;xlabel('east');ylabel('north');zlabel('alt');
        h_ekf_position  = plot3(EP_value(1:index), NP_value(1:index), alt_value(1:index), 'color', 'r','LineWidth', 2);hold on;grid on;
        set(h_ekf_position,'YDataSource','NP_value(1:index)');
        set(h_ekf_position,'XDataSource','EP_value(1:index)');
        set(h_ekf_position,'ZDataSource','alt_value(1:index)');
        set(h_noisy_updates,'YDataSource','ems_data.north_noisy(1:20:index)');
        set(h_noisy_updates,'XDataSource','ems_data.east_noisy(1:20:index)');
        set(h_noisy_updates,'ZDataSource','ems_data.h_noisy(1:20:index)');
        h_vehicle = patch('Faces',faces,'Vertices',CubePoints,'FaceVertexCData',hsv(10),'FaceColor','flat');
        set(h_vehicle,'XData',CubeXData,'YData',CubeYData,'ZData',CubeZData);
        xlabel('EP(m)');ylabel('NP(m)');zlabel('Height(m)');
        %view(-70,30);
        if strcmp(data_type,'sim')
            axis([min(ems_data.east) max(ems_data.east) min(ems_data.north) max(ems_data.north) 650 760]);
        else
            axis([min(ems_data.east) max(ems_data.east) min(ems_data.north) max(ems_data.north) 50 150]);
        end
        legend('Ground Truth','Noisy GPS','EKF');

    else
        if enable_animation == 1
            if (mod(index,sampling_freq) == 0)
                refreshdata(h_ekf_position);
                refreshdata(h_noisy_updates);
                set(h_vehicle,'XData',CubeXData,'YData',CubeYData,'ZData',CubeZData);
                delay(0.05);
            end
        end
    end
    fprintf('processing epoch %d/%d\n',index,num_of_samples);
end

%--> Refreshall figures
refreshdata(h_heading);
refreshdata(h_pitch);
refreshdata(h_roll);
refreshdata(h_Ve);
refreshdata(h_Vn);
refreshdata(h_Vu);
refreshdata(h_east);
refreshdata(h_north);
refreshdata(h_alt);
refreshdata(h_gyro_bias_x);
refreshdata(h_gyro_bias_y);
refreshdata(h_gyro_bias_z);
refreshdata(h_accel_bias_x);
refreshdata(h_accel_bias_y);
refreshdata(h_accel_bias_z);
refreshdata(h_ekf_position);
refreshdata(h_noisy_updates);
set(h_vehicle,'XData',CubeXData,'YData',CubeYData,'ZData',CubeZData);
delay(0.05);

%--> Plot error values
v_east_ref_vector       = ems_data.vel_N(:,1)';
v_north_ref_vector      = ems_data.vel_N(:,2)';
v_up_ref_vector         = ems_data.vel_N(:,3)';
p_east_ref_vector       = ems_data.east';
p_north_ref_vector      = ems_data.north';
alt_ref_vector          = ems_data.h';
roll_ref_vector         = ems_data.roll';
pitch_ref_vector        = ems_data.pitch';
heading_ref_vector      = ems_data.heading';
fprintf('V east error(m/s) = %.10f\n', sqrt(mean((v_east_ref_vector'-ve_value).^2)));
fprintf('V north error(m/s) = %.10f\n', sqrt(mean((v_north_ref_vector'-vn_value).^2)));
fprintf('V up error(m/s) = %.10f\n', sqrt(mean((v_up_ref_vector'-vu_value).^2)));
fprintf('East error(m) = %.10f\n', sqrt(mean((p_east_ref_vector'-EP_value).^2)));
fprintf('North error(m) = %.10f\n', sqrt(mean((p_north_ref_vector'-NP_value).^2)));
fprintf('Alt error(m) = %.10f\n', sqrt(mean((alt_ref_vector'-alt_value).^2)));
fprintf('Roll error(deg) = %.10f\n', sqrt(mean((roll_ref_vector'-Euler_roll_value).^2))*R2D);
fprintf('Pitch error(deg) = %.10f\n', sqrt(mean((pitch_ref_vector'-Euler_pitch_value).^2))*R2D);
fprintf('Heading error(deg) = %.10f\n', sqrt(mean((heading_ref_vector'-Euler_heading_value).^2))*R2D);

figure;title('Orientation Errors (deg)');
subplot(3,1,1);
plot(ems_data.time, (roll_ref_vector' - Euler_roll_value)*R2D,'r');
title('Roll Errors(deg)');
subplot(3,1,2);
plot(ems_data.time, (pitch_ref_vector' - Euler_pitch_value)*R2D,'r');
title('Pitch Errors(deg)');
subplot(3,1,3);
plot(ems_data.time, (heading_ref_vector' - Euler_heading_value)*R2D,'r');
title('Heading Errors(deg)');

figure;title('Velocity Errors');
subplot(3,1,1);
plot(ems_data.time, v_east_ref_vector' - ve_value,'r');
title('Ve Errors(m/s)');
subplot(3,1,2);
plot(ems_data.time, v_north_ref_vector' - vn_value,'r');
title('Vn Errors(m/s)');
subplot(3,1,3);
plot(ems_data.time, v_up_ref_vector' - vu_value,'r');
title('Vu Errors(m/s)');

figure;title('Position Errors');
subplot(3,1,1);
plot(ems_data.time, p_east_ref_vector' - EP_value,'r');
title('East Errors(m)');
subplot(3,1,2);
plot(ems_data.time, p_north_ref_vector' - NP_value,'r');
title('North Errors(m)');
subplot(3,1,3);
plot(ems_data.time, alt_ref_vector' - alt_value,'r');
title('Alt Errors(m)');

figure;
plot(ems_data.time, ems_data.east_noisy, 'r');hold on;grid on;
plot(ems_data.time, ems_data.east, 'b');
plot(ems_data.time, EP_value(1:end), 'k');
xlabel('time(s)');ylabel('east');title('East Position(m)');
legend('Noisy Measurements','Reference','EKF Solution');

figure;
plot(ems_data.time, ems_data.north_noisy, 'r');hold on;grid on;
plot(ems_data.time, ems_data.north, 'b');
plot(ems_data.time, NP_value(1:end), 'k');
xlabel('time(s)');ylabel('north');title('North Position(m)');
legend('Noisy Measurements','Reference','EKF Solution');

figure;
plot(ems_data.time, ems_data.h_noisy, 'r');hold on;grid on;
plot(ems_data.time, ems_data.h, 'b');
plot(ems_data.time, alt_value(1:end), 'k');
if strcmp(data_type,'sim')
    axis([0 500 650 760]);
else
    axis([0 600 50 150]);
end
xlabel('time(s)');ylabel('Alt');title('Alt Position(m)');
legend('Noisy Measurements','Reference','EKF Solution');
