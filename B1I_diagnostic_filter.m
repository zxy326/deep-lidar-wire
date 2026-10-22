classdef B1I_diagnostic_filter < filter_single
% Instrument the existing filter without changing its propagation/update math.
properties
    beforeCorrection=[]
    afterCorrection=[]
    correctionState=[]
end
methods
    function obj=B1I_diagnostic_filter(para)
        obj@filter_single(para);
    end
    function correct(obj,x)
        % INS.solve(...,1) updates quat but leaves its cached Euler angles stale.
        [yaw,pitch,roll]=quat2angle(obj.quat);
        obj.beforeCorrection=[obj.pos,obj.vel,[yaw,pitch,roll]*180/pi];
        correct@filter_single(obj,x);
        obj.afterCorrection=[obj.pos,obj.vel,obj.att];
        obj.correctionState=x;
    end
end
end
