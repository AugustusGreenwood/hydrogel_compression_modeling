function get_lambda_dict()::Dict{Int64, Float64}
    return Dict{Int64, Float64}(
        5 => 14.137166941154069,
        56 => 174.35839227423352,
        35 => 108.38494654884786,
        55 => 171.2167996206437,
        60 => 186.92476288859268,
        30 => 92.6769832808989,
        32 => 98.96016858807849,
        6 => 17.27875959474386,
        67 => 208.91591146372124,
        45 => 139.8008730847458,
        73 => 227.76546738526,
        64 => 199.49113350295187,
        90 => 281.1725424962865,
        4 => 10.995574287564276,
        13 => 39.269908169872416,
        54 => 168.07520696705393,
        63 => 196.34954084936206,
        86 => 268.6061718819273,
        91 => 284.3141351498763,
        62 => 193.20794819577227,
        58 => 180.64157758141312,
        52 => 161.79202165987434,
        12 => 36.12831551628262,
        28 => 86.39379797371932,
        75 => 234.0486526924396,
        23 => 70.68583470577035,
        92 => 287.45572780346606,
        41 => 127.23450247038662,
        43 => 133.5176877775662,
        11 => 32.98672286269283,
        36 => 111.52653920243766,
        68 => 212.05750411731103,
        69 => 215.19909677090084,
        98 => 306.3052837250048,
        82 => 256.0398012675681,
        85 => 265.4645792283375,
        39 => 120.95131716320704,
        84 => 262.32298657474774,
        77 => 240.33183799961918,
        7 => 20.420352248333657,
        25 => 76.96902001294993,
        95 => 296.88050576423547,
        71 => 221.4822820780804,
        66 => 205.77431881013146,
        76 => 237.19024534602937,
        34 => 105.24335389525807,
        50 => 155.50883635269477,
        59 => 183.7831702350029,
        93 => 290.59732045705584,
        2 => 4.71238898038469,
        10 => 29.845130209103033,
        18 => 54.97787143782138,
        26 => 80.11061266653972,
        27 => 83.25220532012952,
        42 => 130.37609512397643,
        87 => 271.7477645355171,
        100 => 312.5884690321844,
        79 => 246.61502330679875,
        16 => 48.69468613064179,
        20 => 61.261056745000964,
        81 => 252.89820861397834,
        19 => 58.119464091411174,
        49 => 152.36724369910496,
        44 => 136.659280431156,
        9 => 26.703537555513243,
        31 => 95.81857593448869,
        74 => 230.9070600388498,
        61 => 190.0663555421825,
        29 => 89.5353906273091,
        94 => 293.7389131106457,
        46 => 142.94246573833559,
        57 => 177.4999849278233,
        70 => 218.34068942449062,
        21 => 64.40264939859075,
        38 => 117.80972450961724,
        88 => 274.8893571891069,
        78 => 243.47343065320896,
        72 => 224.6238747316702,
        24 => 73.82742735936014,
        8 => 23.561944901923447,
        17 => 51.83627878423159,
        37 => 114.66813185602744,
        1 => 1.5707963267948966,
        53 => 164.93361431346415,
        22 => 67.54424205218055,
        47 => 146.08405839192537,
        83 => 259.1813939211579,
        99 => 309.4468763785946,
        89 => 278.0309498426967,
        14 => 42.411500823462205,
        3 => 7.853981633974483,
        80 => 249.75661596038856,
        96 => 300.02209841782525,
        51 => 158.65042900628455,
        33 => 102.10176124166827,
        40 => 124.09290981679683,
        48 => 149.22565104551518,
        15 => 45.553093477052,
        65 => 202.63272615654165,
        97 => 303.16369107141503)
end

function calculate_lambda(n::Int64)::Float64
    return π * (2*n - 1) / 2
end


function get_displacement_summation(t::Float64, sys::System)::Float64
    lambda_dict = get_lambda_dict()
    summation = 0
    for n in 1:100
        λ = lambda_dict[n]
        summation += exp(-λ^2 * t / sys.tg) / λ^2
    end
    return summation
end

function calculate_top_axial_displacement(t::Float64, sys::System)::Float64
    return sys.F * sys.h / (sys.A * sys.Ha) * (1 - 2 * get_displacement_summation(t, sys))
end


function plot_normalized_strain(t::Vector{Float64}, sys::System)
    u = calculate_top_axial_displacement.(t, sys)
    normalized_strain = u * sys.Ha * sys.A / (sys.F * sys.h)

    plot(t/sys.tg, normalized_strain, label="Theory", size=(1200,800))
end


function plot_displacement_validation(file::String, sys::System; plot_options...)
    febio_data = get_febio_data_as_matrix(file)
    febio_time, febio_disp = febio_data[:,1], febio_data[:,2]
    febio_disp_norm = febio_disp / sys.h * sys.Ha * sys.A / sys.F

    plot_normalized_strain(febio_time, sys)    

    time_norm = febio_time / sys.tg

    scatter!(time_norm, febio_disp_norm; plot_options...)
end


function plot_displacement_validation!(file::String, sys::System; plot_options...)
    febio_data = get_febio_data_as_matrix(file)
    febio_time, febio_disp = febio_data[:,1], febio_data[:,2]
    febio_disp_norm = febio_disp / sys.h * sys.Ha * sys.A / sys.F

    time_norm = febio_time / sys.tg

    scatter!(time_norm, febio_disp_norm; plot_options...)
end



function output_displacement_validation(name::String, file::String, sys::System)
    febio_data = get_febio_data_as_matrix(file)
    febio_time, febio_disp = febio_data[:,1], febio_data[:,2]
    febio_disp_norm = febio_disp / sys.h * sys.Ha * sys.A / sys.F

    u = calculate_top_axial_displacement.(febio_time, sys)
    normalized_strain = u * sys.Ha * sys.A / (sys.F * sys.h)

    time_norm = febio_time / sys.tg

    writedlm("$name", [time_norm normalized_strain febio_disp_norm], ',')
end

