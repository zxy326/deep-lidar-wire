%% 北斗 B1I 单天线深组合（独立入口，不修改原 GPS 程序）
clear
clc
root = fileparts(mfilename('fullpath'));
addpath(genpath(root));
cfg = B1I_deep_config();

%% 主机参数与初始条件
Ts = cfg.processSeconds;             % 总处理时长，秒，包含初始化段
To = cfg.startSeconds;               % IQ 文件处理起点，秒
svList = cfg.prns;                   % 北斗 PRN 列表
p0 = [];                            % 可选初始位置 [纬度 经度 高度]；空表示独立定位初始化
psi0 = cfg.deep.initialAttitudeDeg(1);% 初始航向，度；未知时为 NaN，需实测填写
pitch0 = cfg.deep.initialAttitudeDeg(2);
roll0 = cfg.deep.initialAttitudeDeg(3);
arm = cfg.deep.armBodyMeters;        % IMU 指向天线的杆臂，前右下机体系，米
gyro0 = cfg.deep.gyroReferenceDegps;  % 运动检测参考，不将未知运动的均值当零偏

%% 时间同步与 IMU 安装信息
% 在 B1I_deep_config.m 中填写，不能随意设为已同步：
% cfg.sync.confirmed、imuReferenceSeconds、imuToIqOffsetSeconds、imuClockScale
% cfg.sync.imuAxesToBody：列向量 body = R * sensor
% IQ秒数 = (IMU时间戳 - imuReferenceSeconds) * imuClockScale + 偏移
% cfg.imu.specificForceConfirmed：确认加速度为比力，未去除重力响应
% cfg.imu.topic：bag 内明确的 sensor_msgs/Imu 主题

%% 选择 IMU 数据文件
[file, path] = uigetfile({'*.bag;*.mat;*.csv','IMU 数据 (*.bag, *.mat, *.csv)'}, ...
    '选择 IMU 数据文件',cfg.imuFile);
if isequal(file,0), error('B1I:Cancelled','未选择 IMU 文件。'); end
cfg.imuFile = fullfile(path,file);
[~,~,extension] = fileparts(file);
switch lower(extension)
    case '.bag', cfg.imu.format = 'ros1bag';
    case '.mat', cfg.imu.format = 'mat';
    case '.csv', cfg.imu.format = 'csv';
end
% MAT 变量 imu / 无表头 CSV：时间、wx wy wz、ax ay az，共七列。
imu = BDS.B1I.load_imu(cfg);          % 输出已同步的 IQ 秒数、rad/s、m/s^2
imuN = size(imu,1);

%% 选择北斗 B1I 中频 IQ 数据文件
[file, path] = uigetfile({'*.dat;*.bin','IQ 数据 (*.dat, *.bin)'}, ...
    '选择北斗 B1I 中频数据文件',cfg.iqFile);
if isequal(file,0), error('B1I:Cancelled','未选择 IQ 文件。'); end
data_file = fullfile(path,file);
sampleFreq = cfg.sampleFreq;         % 当前录制：5 MHz，需与采集配置一致
sampleOffset = round(To*sampleFreq); % 复数采样点偏移
blockSize = round(sampleFreq*0.001);  % 名义 1 ms；实际通道按码周期读取
msToProcess = ceil(Ts*1000);
% 当前 IQ 为 float32 I/Q 交织，每复数采样 8 字节，不是 GPS 的 int16。
% 文件由 BB1I_S / DeepChannel 打开和关闭，不重复创建外部 fileID。

%% 接收机配置
receiver_conf = cfg;
receiver_conf.iqFile = data_file;
receiver_conf.startSeconds = To;
receiver_conf.processSeconds = Ts;
receiver_conf.sampleFreq = sampleFreq;
receiver_conf.prns = svList;
receiver_conf.deep.initialAttitudeDeg = [psi0,pitch0,roll0];
receiver_conf.deep.armBodyMeters = arm;
if ~isempty(p0), receiver_conf.initialEcef = lla2ecef(p0); end
BDS.B1I.validate_deep_setup(receiver_conf,imu);
% 不从文件名猜 GPS 时间：B1I 电文提供 BDT，IMU 使用明确的时间映射。
% 不下载 GPS 历书；B1I 初始化过程中自行解码 D1/D2 星历。

%% 导航滤波器参数（与原 GPS 主脚本对应）
receiver_conf.deep.gyroReferenceDegps = gyro0;
% 与非交互入口共用唯一参数构造，包含物理伪距率模型开关。
para = BDS.B1I.deep_parameters(receiver_conf,[0,0,0],[0,0,0]);

%% 创建北斗 B1I 接收机对象
nCoV = BB1I_S(receiver_conf);
f = [];                            % 进度条句柄
tic
try
    %% 标量捕获、跟踪、电文和初始定位
    % 完成后 state=1；这是离线初始化，不是实时硬件接收机。
    nCoV.initialize();
    initialTime = nCoV.initialTime;
    u0 = BDS.B1I.imu_at(imu,initialTime,receiver_conf);

    %% 进入深组合（与 GPS 的 state==1 分支对应）
    if nCoV.state==1
        nCoV.imu_input(initialTime,u0); % 已转换为 rad/s，不重复转单位
        para.p0 = nCoV.pos;
        para.v0 = nCoV.vel;
        nCoV.navFilter = B1I_navigation_filter(para);
        nCoV.vectorMode = 2;
        nCoV.channel_vector();        % 初始化杆臂/钟差，切换到 state=3
    end

    %% 创建进度条
    f = waitbar(0,'北斗 B1I + IMU 深组合');
    replayMs = ceil((nCoV.endTime-initialTime)*1000);

    %% 接收机运行：1 ms 调度，按 IMU 周期融合
    for t=1:replayMs
        targetTime = min(initialTime+t*0.001,nCoV.endTime);
        % 与非交互入口共用调度；严格先处理到融合时刻并反馈，再继续 IQ。
        BDS.B1I.advance_receiver(nCoV,targetTime,imu);
        if mod(t,1000)==0 && isgraphics(f)
            waitbar(t/replayMs,f,sprintf('IQ %.2f / %.2f 秒',nCoV.ta,nCoV.endTime));
        end
        if nCoV.state==0, break; end
    end
    if nCoV.state~=0, nCoV.stop('Requested interval completed'); end

    %% 整理结果
    nCoV.clean_storage();
    b1iDeepResult = nCoV.get_result();
    toc
    if isgraphics(f), close(f); end

    %% 保存 B1I 星历和深组合结果（每次新建目录，不覆盖）
    if ~isfolder(receiver_conf.outputDir), mkdir(receiver_conf.outputDir); end
    runDir = tempname(receiver_conf.outputDir); mkdir(runDir);
    ephemeris_file = fullfile(runDir,'b1i_ephemeris.mat');
    nCoV.save_ephemeris(ephemeris_file);
    b1iDeepResult.runDir = runDir;
    result = b1iDeepResult;          % 两个入口统一保存变量；工作区保留原名称。
    save(fullfile(runDir,'deep_result.mat'),'result','receiver_conf','para','data_file','-v7.3');

    %% 画融合轨迹、速度和量测更新情况
    if receiver_conf.showPlot, nCoV.plot_result(); end
    fprintf('处理结束：%s\n结果目录：%s\n',b1iDeepResult.stopReason,runDir);
catch err
    if ~isempty(f) && isgraphics(f), close(f); end
    delete(nCoV);                   % 异常也关闭各通道 IQ 文件
    rethrow(err);
end
delete(nCoV);                       % 结果保留在 b1iDeepResult，不保存活动文件句柄
