%   Delay module to refresh plots

%   Symbolic definistions of system dynamics, error models, and EKF models

%   This is a basic INS/GNSS loosely-coupled system using EKF
%   This is only a DEMO with basic updates (position/velocity) are applied
%   More advanced updates such as nonholonomic constraints or zero-speed
%   are NOT implemented in this DEMO. The purpose of the DEMO is to
%   demonstrate the basics of EKF in a basic INS/GNSS fusion scenario

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

function delay(T)
if(T > 0) %end condition
    t = timer('timerfcn','delay(0)','StartDelay',T);
    start(t);
    wait(t);
end