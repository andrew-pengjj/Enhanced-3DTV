if strcmp('on', get(hpan, 'Enable'))
    set(hpan, 'Enable', 'off');
else
    set(hpan, 'Enable', 'on');
    set(hzoom, 'Enable', 'off');
    set(bzoom, 'Value', get(bzoom, 'Min'));
end