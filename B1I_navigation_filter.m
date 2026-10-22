classdef B1I_navigation_filter < B1I_diagnostic_filter
% B1I-local version of the original filter run(), with physical range-rate
% scale. Common/GPS classes are untouched. Disabled mode calls the original.
properties
 physicalRangeRate=false
end
methods
 function obj=B1I_navigation_filter(para)
  obj@B1I_diagnostic_filter(para);
  if isfield(para,'physicalRangeRate'), obj.physicalRangeRate=para.physicalRangeRate; end
 end
        function run(obj, imu, sv, indexP, indexV)
            if nargin==2
                flag = 0;
            else
                flag = 1;
            end
            rateOnly=flag==1 && ~any(indexP) && any(indexV);
            if ~obj.physicalRangeRate && ~rateOnly
                if nargin==2, run@filter_single(obj,imu); else, run@filter_single(obj,imu,sv,indexP,indexV); end
                return;
            end
            r2d = 180/pi;
            c = 299792458;
            dt = obj.T;
            wbo = imu(1:3);
            obj.motion.run(wbo*r2d, obj.vel);
            obj.wdot = obj.wdotCal.run(wbo);
            wbd = obj.wbDelay.push(wbo);
            wbd = wbd - obj.bias(1:3);
            imu = imu - obj.bias;
            wb1 = imu(1:3);
            fb1 = imu(4:6);
            obj.solve(imu, 1);
            obj.dtr = obj.dtr + obj.dtv*dt;
            Cnb = quat2dcm(obj.quat);
            Cbn = Cnb';
            fn = fb1*Cnb;
            winn = obj.geogInfo.wien + obj.geogInfo.wenn;
            A = zeros(17);
            A(1:3,12:14) = -Cbn;
            A(4:6,1:3) = antisym(fn);
            A(4:6,15:17) = Cbn;
            A(7,4) = 1;
            A(8,5) = 1;
            A(9,6) = 1;
            A(10,11) = 1;
            Phi = eye(17) + A*dt + (A*dt)^2/2;
            P1 = Phi*obj.P*Phi' + obj.Q;
            X = zeros(17,1);
            measureFlag = 0;
            if flag==1 && (any(indexP)||any(indexV)) && obj.accJump.state==0
                measureFlag = 1;
                n1 = sum(indexP);
                n2 = sum(indexV);
                rs = sv(:,1:3);
                vs = sv(:,4:6);
                rho = sv(:,7);
                rhodot = sv(:,8);
                R_rho = sv(:,9);
                R_rhodot= sv(:,10);
                if obj.delay==0
                    [rho0, rhodot0, rspu, Cen] = rho_rhodot_cal_geog(rs, vs, obj.pos, obj.vel);
                else
                    [pos, vel] = delayComp(obj.pos, obj.vel, obj.acc, obj.geogInfo.Cn2g, obj.delay);
                    [rho0, rhodot0, rspu, Cen] = rho_rhodot_cal_geog(rs, vs, pos, vel);
                end
    
                if obj.physicalRangeRate
                    cm = BDS.B1I.rate_scale(rs,vs,lla2ecef(obj.pos));
                else
                    cm = 1-sum(rspu.*vs,2)/c;
                end
                En = rspu*Cen';
                rho = rho - obj.dtr*c;
                rhodot = rhodot - obj.dtv*c - obj.windupFlag*wb1(3)*0.030286178664972;
                ran = Cbn*obj.arm';
                rho = rho - En*ran;
                van = Cbn*cross(wb1,obj.arm)';
                rhodot = rhodot - En*van;
                E = En(indexP,:);
                H1 = zeros(n1,17);
                H1(:,7:9) = E;
                H1(:,10) = -1;
                Z1 = rho0(indexP) - rho(indexP);
                R1 = diag(R_rho(indexP));
                H2 = []; Z2 = []; R2 = [];
                if n2>0
                    E = En(indexV,:);
                    H2 = zeros(n2,17);
                    H2(:,4:6) = E;
                    H2(:,11) = -cm(indexV);
                    Z2 = rhodot0(indexV) - rhodot(indexV).*cm(indexV);
                    if obj.motion.state==0
                        R2 = diag(R_rhodot(indexV));
                    else
                        R2 = diag(R_rhodot(indexV)*1);
                    end
                end
                H3 = []; Z3 = []; R3 = [];
                if obj.motion.state==0
                    H3 = zeros(3,17);
                    H3(:,12:14) = eye(3);
                    Z3 = (wbd-winn*Cbn)';
                    R3 = diag([1;1;1]*obj.Rwb);
                end
                H = [H1; H2; H3];
                Z = [Z1; Z2; Z3];
                R = blkdiag(R1, R2, R3);
                K = P1*H' / (H*P1*H'+R);
                X = K*Z;
                R_diag = diag(R);
                R_sqrt_diag = sqrt(R_diag);
                gamma = 1.3;
                for k=1:2
                    Psi_R_diag = HuberWeight((Z-H*X)./R_sqrt_diag, gamma);
                    R0 = diag(R_diag./Psi_R_diag);
                    K = P1*H' / (H*P1*H'+R0);
                    X = K*Z;
                end
                P1 = (eye(17)-K*H)*P1;
            elseif obj.motion.state==0
                measureFlag = 2;
                H = zeros(3,17);
                H(:,12:14) = eye(3);
                Z = (wbd-winn*Cbn)';
                R = diag([1;1;1]*obj.Rwb);
                K = P1*H' / (H*P1*H'+R);
                X = K*Z;
                P1 = (eye(17)-K*H)*P1;
            end
            if measureFlag~=0
                Y = [];
                if obj.motion.state==0
                    Ysub = zeros(1,17);
                    Ysub(1,1) = Cnb(1,1)*Cnb(1,3);
                    Ysub(1,2) = Cnb(1,2)*Cnb(1,3);
                    Ysub(1,3) = -(Cnb(1,1)^2+Cnb(1,2)^2);
                    Y = [Y; Ysub];
                    Ysub = zeros(2,17);
                    Ysub(1,15) = 1;
                    Ysub(2,16) = 1;
                    Y = [Y; Ysub];
                end
                if ~isempty(Y)
                    X = X - P1*Y'/(Y*P1*Y')*Y*X;
                end
            end
            obj.P = (P1+P1')/2;
            X = X';
            obj.correct(X(1:9));
            obj.dtr = obj.dtr + X(10)/c;
            obj.dtv = obj.dtv + X(11)/c;
            obj.bias(1:3) = obj.bias(1:3) + X(12:14);
            obj.bias(4:6) = obj.bias(4:6) + X(15:17);
        end
end
end
