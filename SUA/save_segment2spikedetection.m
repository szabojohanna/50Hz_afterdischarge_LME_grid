function save_segment2spikedetection(cnt_path,cnt_filename,time_start_sec, time_end_sec,chan);


    nsT = nswiew;
    hT = guidata(nsT);
    
    
    hT = nswiew('inport_menu_Callback',hT, [], hT, { [cnt_path '\'], cnt_filename});
    
    guidata(nsT, hT);
    hT =guidata(nsT);

    time2load_sec = diff([time_start_sec time_end_sec]);
    hT = nswiew('inport',time_start_sec,time2load_sec,hT);
    guidata(nsT,hT)
    hT = guidata(nsT);
        
    data = hT.data(:,chan);
    sr = hT.srate;
    [savedir, fnm2] = save_fnm(cnt_path,cnt_filename,time_start_sec, time_end_sec,chan);
    
    try
    save(fullfile(savedir,[fnm2 '.mat']),'data','sr');
    catch
        save(fullfile(savedir,[fnm2 '.mat']),'data','sr','-v7.3');
    end
close(nsT);
end