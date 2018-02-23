%   convert from DCM to quaternion

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

function quat_vector = convert_dcm_to_quat( C )
    Tr = trace(C);
    Pa = 1 + Tr;
    Pb = 1 + 2*C(1,1)-Tr;
    Pc = 1 + 2*C(2,2)-Tr;
    Pd = 1 + 2*C(3,3)-Tr;

    [m,i] = max([Pa Pb Pc Pd]);
    
    switch(i)
        case 1 
        a = 0.5*sqrt(Pa);
        b = ( C(3,2)-C(2,3) )/(4*a);
        c = (C(1,3) - C(3,1))/(4*a); 
        d = (C(2,1) - C(1,2))/(4*a);
        case 2 
        b = 0.5*sqrt(Pb);
        c = ( C(2,1)+C(1,2) )/(4*b); 
        d = (C(1,3) + C(3,1))/(4*b); 
        a = (C(3,2) - C(2,3))/(4*b);
        case 3 
        c = 0.5*sqrt(Pc);
        d = ( C(3,2)+C(2,3) )/(4*c);
        a = (C(1,3) - C(3,1))/(4*c); 
        b = (C(2,1) + C(1,2))/(4*c);
        case 4 
        d = 0.5*sqrt(Pd);
        a = ( C(2,1)-C(1,2) )/(4*d); 
        b = (C(1,3) + C(3,1))/(4*d); 
        c = (C(3,2) + C(2,3))/(4*d);
    end
    
    if a <= 0
        a = -a;
        b = -b;
        c = -c;
        d = -d;
    end

    quat_vector = [a;b;c;d];
    
end