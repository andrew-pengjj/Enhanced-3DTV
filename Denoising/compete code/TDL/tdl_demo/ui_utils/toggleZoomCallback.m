if strcmp('on', get(hzoom, 'Enable'))
    set(hzoom, 'Enable', 'off');
else
    set(hzoom, 'Enable', 'on');
    set(hpan, 'Enable', 'off');
    set(bpan, 'Value', get(bpan, 'Min'));
end