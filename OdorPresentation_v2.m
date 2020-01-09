function OdorPresentation_v2(COMPort, bee, sess)

%-------------------------------------------------------

% For example:

%   OdorPresentation_v2('COM4', 400, 2)



global ex





% gui initialization

init_plugin_gui()

init_RunPauseButton();

init_StopStartButton();





ex.tline = [];

ex.COMPort = COMPort;

ex.bee = bee;

ex.sess = sess;

ex.protocol = 'OdorPresentation';



%ex.ard = s;

tic



end



function start_session()

global ex tr



ex.fpath = fullfile('C:', 'Experiments', sprintf('bee_%03d',ex.bee));

ex.fname = sprintf('m%03ds%02d.txt',ex.bee, ex.sess);

ex.mname = sprintf('m%03ds%02d.mat',ex.bee, ex.sess);



if ~isdir(ex.fpath)

    mkdir(ex.fpath)

end



fn  = fullfile(ex.fpath, ex.fname);



fid = fopen(fn, 'wt');

fprintf(fid, '# Recording date : %s \n', datestr(now, 'dd-mmm-yyyy'));

fprintf(fid, '# Matlab Function name - %s \n', mfilename);

fprintf(fid, '# Protocol - %s \n', ex.protocol);

fprintf(fid, '#\n');

fprintf(fid, '# start : %s \n\n', datestr(now, 'HH:MM:SS'));

fprintf(fid, 'tr\t t0\t Odor\t OdorName\t OdorConc\t ITI\n');

fprintf(fid, '# -----------------------------------------------------------------------------\n');



tr = [];

ex.fid = fid;

ex.tr_num = 0;

ex.ind = 0;

ex.rand_arr = randperm(500);



% setting up serial port

ex.ard = CreateArduinoComm(ex.COMPort);



end



function end_session()

global ex

% finishing behavior file

fid = ex.fid;

fprintf(fid, '# -----------------------------------------------------------------------------\n');

fprintf(fid, '# stop : %s \n', datestr(now, 'HH:MM:SS'));

fclose (fid);





% disconnecting from serial port

ReleaseArduino(ex.ard)

disp('Serial communication is over')

ex.tr_num = 0;



end



function s = CreateArduinoComm(COMPort)

 

delete(instrfind) 

s=serial(COMPort);

set(s, 'Timeout', 30); %Sets the timeout to 30 seconds

 

 

s.baudrate=115200;

s.flowcontrol='none';

s.inputbuffersize=100;

s.bytesavailablefcnmode = 'terminator';

 

 

set(s,'Terminator','CR/LF'); %/LF

set(s,'DataBits',8);

set(s,'StopBits',2);

 

 

fopen(s);

pause(0.1);



disp('Serial communication is ready')

 

end

 

function ReleaseArduino(s)

 

fclose(s);

delete(s)

clear s

 

delete(instrfind) %to clear stack overflow

 

end



function init_plugin_gui()

% gui parameter initialization

% main olfactometer figure

global gp

 

    % color assignement for odor, water and sound

    gp.fig_color    = [0.8, 0.9, 0.8];



    gp.scr_size         = get(0,'ScreenSize'); % [left, bottom, width, height]:

    gp.fig_size         = [200, 200];

    gp.fig_pos          = [200, 2*gp.scr_size(3)/4-200, gp.fig_size];

    gp.axes_pos         = [0.4, 0.06, 0.54, 0.88 ];

    gp.but_x_pos        = [70 190 190 100];

    gp.but_y_pos        = [100 50 330 305];

    gp.button_size      = [80 30];

    gp.radio_size       = [60 30];

    gp.button_color     = [0.5 0.7 0.9];

    gp.background_color = [0.84,0.91,0.85];



    gp.fig  = figure(20);

    clf;

    set(gp.fig, 'Position',    gp.fig_pos, ...

                'Resize',      'off', ...

                'NumberTitle', 'off', ...

                'Name',        'Detect_v1', ...

                'Color',       'gp.fig_color');



    drawnow

end



function init_RunPauseButton()

global gp

pos = [gp.but_x_pos(1) gp.but_y_pos(1)];

gp.pause = uicontrol(gp.fig, 'Style',              'ToggleButton', ...

                             'position',           [pos, gp.button_size], ...

                             'BackgroundColor',    gp.button_color,...

                             'string',             'Run',...

                             'callback',           @pauseButton);

    

    function pauseButton(~, ~)

    

    pause_flag=get(gp.pause,'Value');

    start_flag = get(gp.stop,'Value');

%     fprintf('Is Start button pressed - %1d\n',start_flag)

    if (pause_flag == 1)&&(start_flag)

        set(gp.pause,'BackgroundColor',[0.4 1 0.0]);

        set(gp.pause,'Value', 1);

        set(gp.pause,'string', 'Pause');

        protocol();

    else

        set(gp.pause,'BackgroundColor', gp.button_color);

        set(gp.pause,'string', 'Run');

        set(gp.pause,'Value', 0);

    end



    end

end



function init_StopStartButton()

global gp

pos = [gp.but_x_pos(1) gp.but_y_pos(2)];

gp.stop = uicontrol(gp.fig, 'Style',              'ToggleButton', ...

                            'position',           [pos, gp.button_size], ...

                            'BackgroundColor',    gp.button_color,...

                            'string',             'Start session',...

                            'callback',           @stopButton);

    

    function stopButton(~, ~)

    

    start_flag=get(gp.stop,'Value');

    run_flag = get(gp.pause,'Value');

    if (start_flag == 1)&&(run_flag==0)

        set(gp.stop,'BackgroundColor',[0.4 1 0.0]);

        set(gp.stop,'string', 'Stop session');

        set(gp.stop,'Value', 1);

        start_session();

    elseif (start_flag==0)&&(run_flag==0)

        set(gp.stop,'BackgroundColor', gp.button_color);

        set(gp.stop,'string', 'Start session');

        set(gp.stop,'Value', 0);

        end_session();

        

    elseif (start_flag==0)&&(run_flag==1)

        set(gp.stop,'BackgroundColor', [0.4 1 0.0]);

        set(gp.stop,'string', 'Stop session');

        set(gp.stop,'Value', 1);

    end



    end

end



function protocol()

global ex gp tr

ex.tline = [];



% setting up session parameters

ITI       = uint16(7000);   % interval between trials in msec

Ntrials   = 10;   % per each stimulus

odorDur   = 1000; % in msec



ss(1) = struct('odorNum', 5 , 'odorName', 'Linalool', 'odorConc', 100);

ss(2) = struct('odorNum', 5 , 'odorName', 'Linalool', 'odorConc',  33);

ss(3) = struct('odorNum', 5 , 'odorName', 'Linalool', 'odorConc',  11);

ss(4) = struct('odorNum', 6 , 'odorName', 'BP', 'odorConc', 100);

ss(5) = struct('odorNum', 6 , 'odorName', 'BP', 'odorConc',  33);

ss(6) = struct('odorNum', 6 , 'odorName', 'BP', 'odorConc',  11);

ss(7) = struct('odorNum', 7 , 'odorName', 'Ethanol', 'odorConc',  100);



% do not modify below



ss = repmat(ss, 1, Ntrials);

ex.n_tr = length(ss);

randVec = randperm(ex.n_tr);



tline = ex.tline;

s = ex.ard;



% setting up session parameters

ex.params.date = datestr(now, 'dd-mmm-yyyy');

ex.params.behavPar = ex.protocol;

ex.params.total_flow = 1000; % in msec

ex.params.N2_flow = 50; % in msec

ex.params.ITI_dur = ITI; % in msec



pause_flag = get(gp.pause,'Value'); % if pause button is pressed



while pause_flag && (ex.tr_num < ex.n_tr)



    ex.ind = ex.ind+1;

    ex.tr_num = ex.tr_num + 1;

    if ex.tr_num==1

        tic

    end

    ex.t0 = round(1e3*toc);

    

    tr(ex.tr_num).date = ex.params.date;

    tr(ex.tr_num).behavPar = ex.params.behavPar;

    tr(ex.tr_num).total_flow = ex.params.total_flow;

    tr(ex.tr_num).N2_flow = ex.params.N2_flow;

    tr(ex.tr_num).ITI_dur = ex.params.ITI_dur;

    

    

    %rand_ind = ceil(length(odorNum)*rand(1));

    

    tr(ex.tr_num).odorName  = ss(randVec(ex.ind)).odorName;

    tr(ex.tr_num).odorNum   = ss(randVec(ex.ind)).odorNum;

    tr(ex.tr_num).odorConc  = ss(randVec(ex.ind)).odorConc;

    tr(ex.tr_num).odorDur   = odorDur;

    tr(ex.tr_num).t0        = ex.t0;

    tr(ex.tr_num).tr_num    = ex.tr_num;

    

    %display in command window

    fprintf('%d : %d\t %d \r', ex.tr_num, tr(ex.tr_num).t0, tr(ex.tr_num).odorNum, tr(ex.tr_num).odorConc);

    

    %sentence that will be sended to Arduino

    s1 = sprintf('%s %d %d %d %d\r', 'params',...

                  tr(ex.tr_num).ITI_dur, tr(ex.tr_num).odorNum, tr(ex.tr_num).odorConc, tr(ex.tr_num).odorDur);

    % sending parameters to Arduino

    fprintf(ex.ard, s1);

 

    while ~strncmp('Done...',tline, 7) && ~strncmp('>Done...',tline, 8)

        if s.BytesAvailable

            %sprintf('There are %d bytes\n', s.BytesAvailable)

            tline = fgetl(s);

            disp(tline)

        end

    end

 

    % saving trial parameters to file

    fprintf(ex.fid, '%d : \t %d\t %4d\t %s\t\t %d\t\t %d \n\n',...

            tr(ex.tr_num).tr_num, tr(ex.tr_num).t0, tr(ex.tr_num).odorNum ,...

            tr(ex.tr_num).odorName, tr(ex.tr_num).odorConc, tr(ex.tr_num).ITI_dur);

    

    save(fullfile(ex.fpath,ex.mname), 'tr');

                

    pause(0.1)

    tline = '';

    pause_flag = get(gp.pause, 'Value');

    

end % end of while



end

