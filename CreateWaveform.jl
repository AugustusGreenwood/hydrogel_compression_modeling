function __create_cycle(t0, y0, A, ramp, dwell)
    return [t0+ramp y0+A; 
            t0+ramp+dwell y0+A;
            t0+2ramp+dwell y0;
            t0+2*(ramp+dwell) y0]
end

function __waveform_center(cycles, t0, y0, A, ramp, dwell)
    waveform = __create_cycle(t0, y0, A, ramp, dwell)
    for i in 1:cycles
        waveform = [waveform; __create_cycle(t0+2i*(ramp+dwell), y0, A, ramp, dwell)]
    end
    return waveform
end

function __waveform_start(A1)
    return [0.0 0.0; 0.5 A1; 1.5 A1]
end


function __waveform_end(t0, A1)
    return [t0+1.0 A1; t0+1.5 0.0]
end

function create_waveform(A1, A2, ramp, dwell, cycles)
    waveform = __waveform_start(A1)
    waveform = [waveform; __waveform_center(cycles, 1.5, A1, A2, ramp, dwell)]
    waveform = [waveform; waveform[end,1]+0.5 0.0]
    return waveform
end




