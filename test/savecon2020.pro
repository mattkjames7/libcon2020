pro savexyz

    x = double(randomu(234827,100)*200 - 100)
    y = double(randomu(2264827,100)*200 - 100)
    z = double(randomu(2264827,100)*50 - 25)

    b = con2020_model_xyz("analytic",x,y,z)
    bx = b[*,0]
    by = b[*,1]
    bz = b[*,2]

    openw, lun, "con2020-analytic-xyz.bin", /get_lun
    writeu, lun, long(100)
    writeu, lun, double(bx)
    writeu, lun, long(100)
    writeu, lun, double(by)
    writeu, lun, long(100)
    writeu, lun, double(bz)
    free_lun, lun

    b = con2020_model_xyz("integral",x,y,z)
    bx = b[*,0]
    by = b[*,1]
    bz = b[*,2]

    openw, lun, "con2020-integral-xyz.bin", /get_lun
    writeu, lun, long(100)
    writeu, lun, double(bx)
    writeu, lun, long(100)
    writeu, lun, double(by)
    writeu, lun, long(100)
    writeu, lun, double(bz)
    free_lun, lun

    b = con2020_model_xyz("hybrid",x,y,z)
    bx = b[*,0]
    by = b[*,1]
    bz = b[*,2]

    openw, lun, "con2020-hybrid-xyz.bin", /get_lun
    writeu, lun, long(100)
    writeu, lun, double(bx)
    writeu, lun, long(100)
    writeu, lun, double(by)
    writeu, lun, long(100)
    writeu, lun, double(bz)
    free_lun, lun

end

pro savertp

    r = double(randomu(234827,100)*100)
    t = double(randomu(2264827,100)*!DPI)
    p = double(randomu(2264827,100)*2*!DPI)

    b = con2020_model_rtp("analytic",r,t,p)
    br = b[*,0]
    bt = b[*,1]
    bp = b[*,2]

    openw, lun, "con2020-analytic-rtp.bin", /get_lun
    writeu, lun, long(100)
    writeu, lun, double(br)
    writeu, lun, long(100)
    writeu, lun, double(bt)
    writeu, lun, long(100)
    writeu, lun, double(bp)
    free_lun, lun

    b = con2020_model_rtp("integral",r,t,p)
    br = b[*,0]
    bt = b[*,1]
    bp = b[*,2]

    openw, lun, "con2020-integral-rtp.bin", /get_lun
    writeu, lun, long(100)
    writeu, lun, double(br)
    writeu, lun, long(100)
    writeu, lun, double(bt)
    writeu, lun, long(100)
    writeu, lun, double(bp)
    free_lun, lun

    b = con2020_model_rtp("hybrid",r,t,p)
    br = b[*,0]
    bt = b[*,1]
    bp = b[*,2]

    openw, lun, "con2020-hybrid-rtp.bin", /get_lun
    writeu, lun, long(100)
    writeu, lun, double(br)
    writeu, lun, long(100)
    writeu, lun, double(bt)
    writeu, lun, long(100)
    writeu, lun, double(bp)
    free_lun, lun
end