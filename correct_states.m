%   Correct States

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

%--> Use the ekf state_correction_vector to correct all states 

lat_value(index+1) = (lat_value(index+1)*D2R + state_correction_vector(2)/(Rm_value(index+1)+alt_value(index+1)))*R2D;
lon_value(index+1) = (lon_value(index+1)*D2R + state_correction_vector(1)/((Rm_value(index+1)+alt_value(index+1))*cos(lat_value(index+1)*D2R)))*R2D;
C_EN_value = C_EN_fun(lat_value(index+1)*D2R , lon_value(index+1)*D2R);
pos_quat_vector = convert_dcm_to_quat (C_EN_value);
a_pos_value(index+1) = pos_quat_vector(1);
b_pos_value(index+1) = pos_quat_vector(2);
c_pos_value(index+1) = pos_quat_vector(3);
d_pos_value(index+1) = pos_quat_vector(4);
normalization_factor = sqrt(a_pos_value(index+1)^2 + b_pos_value(index+1)^2 + c_pos_value(index+1)^2 + d_pos_value(index+1)^2);
a_pos_value(index+1) = a_pos_value(index+1)/normalization_factor;
b_pos_value(index+1) = b_pos_value(index+1)/normalization_factor;
c_pos_value(index+1) = c_pos_value(index+1)/normalization_factor;
d_pos_value(index+1) = d_pos_value(index+1)/normalization_factor;
alt_value(index+1) = alt_value(index+1) + state_correction_vector(3);
ve_value(index+1) = ve_value(index+1) + state_correction_vector(4);
vn_value(index+1) = vn_value(index+1) + state_correction_vector(5);
vu_value(index+1) = vu_value(index+1) + state_correction_vector(6);
b_value(index+1) = b_value(index+1) + state_correction_vector(7);
c_value(index+1) = c_value(index+1) + state_correction_vector(8);
d_value(index+1) = d_value(index+1) + state_correction_vector(9);
a_value(index+1) = a_value(index+1) - ([b_value(index) c_value(index) d_value(index)]*state_correction_vector(7:9))/a_value(index);
gyro_bias_x_value(index+1) = gyro_bias_x_value(index+1) + state_correction_vector(10);
gyro_bias_y_value(index+1) = gyro_bias_y_value(index+1) + state_correction_vector(11);
gyro_bias_z_value(index+1) = gyro_bias_z_value(index+1) + state_correction_vector(12);
acc_bias_x_value(index+1) = acc_bias_x_value(index+1) + state_correction_vector(13);
acc_bias_y_value(index+1) = acc_bias_y_value(index+1) + state_correction_vector(14);
acc_bias_z_value(index+1) = acc_bias_z_value(index+1) + state_correction_vector(15);
gyro_sf_x_value(index+1) = gyro_sf_x_value(index+1) + state_correction_vector(16);
gyro_sf_y_value(index+1) = gyro_sf_y_value(index+1) + state_correction_vector(17);
gyro_sf_z_value(index+1) = gyro_sf_z_value(index+1) + state_correction_vector(18);
acc_sf_x_value(index+1) = acc_sf_x_value(index+1) + state_correction_vector(19);
acc_sf_y_value(index+1) = acc_sf_y_value(index+1) + state_correction_vector(20);
acc_sf_z_value(index+1) = acc_sf_z_value(index+1) + state_correction_vector(21);

% Normalize the quaternions
normalization_factor = sqrt(a_value(index+1)^2 + b_value(index+1)^2 + c_value(index+1)^2 + d_value(index+1)^2);
a_value(index+1) = a_value(index+1)/normalization_factor;
b_value(index+1) = b_value(index+1)/normalization_factor;
c_value(index+1) = c_value(index+1)/normalization_factor;
d_value(index+1) = d_value(index+1)/normalization_factor;
normalization_factor = sqrt(a_pos_value(index+1)^2 + b_pos_value(index+1)^2 + c_pos_value(index+1)^2 + d_pos_value(index+1)^2);
a_pos_value(index+1) = a_pos_value(index+1)/normalization_factor;
b_pos_value(index+1) = b_pos_value(index+1)/normalization_factor;
c_pos_value(index+1) = c_pos_value(index+1)/normalization_factor;
d_pos_value(index+1) = d_pos_value(index+1)/normalization_factor;

% Calculate DCM from Quaternion
C_LB_from_quat_value = C_LB_from_quat_fun(...
    b_value(index+1),...
    c_value(index+1),...
    d_value(index+1));

% Calculate Euler angles from quaternion
Euler_pitch_value(index+1)     = atan2(-C_LB_from_quat_value(3,1),(sqrt(C_LB_from_quat_value(3,2)^2 + C_LB_from_quat_value(3,3)^2)));
Euler_roll_value(index+1)      = atan2(C_LB_from_quat_value(3,2),C_LB_from_quat_value(3,3));
Euler_heading_value(index+1)   = atan2(C_LB_from_quat_value(2,1),C_LB_from_quat_value(1,1));
Rn_value(index+1)  = earth_a/sqrt(1-earth_e2*sin(lat_value(index+1)*D2R)*sin(lat_value(index+1)*D2R));
Rm_value(index+1)  = earth_a*(1-earth_e2)/((1-earth_e2*sin(lat_value(index+1)*D2R)*sin(lat_value(index+1)*D2R))^(1.5));
EP_value(index+1) = (lon_value(index+1)-lon_value(1))*D2R*(Rn_value(index+1)+alt_value(index+1))*cos(lat_value(index+1)*D2R);
NP_value(index+1) = (lat_value(index+1)-lat_value(1))*D2R*(Rm_value(index+1)+alt_value(index+1));