

def tskypy(gb,gl):
    """ Calculate tsky from Haslam table, scale to survey frequency"""
    tskypath = "/Users/aemcewen/GBNCC_paper/jail/tsky.ascii"
    tskylist = []
    with open(tskypath) as f:
        for line in f:
            str_idx = 0
            while str_idx < len(line):
                # each temperature occupies space of 5 chars
                temp_string = line[str_idx:str_idx+5]
                try:
                    tskylist.append(float(temp_string))
                except:
                    pass
                str_idx += 5

    # ensure l is in range 0 -> 360
    b = gb
    if gl < 0.:
        l = 360 + gl
    else:
        l = gl

    # convert from l and b to list indices
    j = b + 90.5
    if j > 179:
        j = 179
    nl = l - 0.5
    if l < 0.5:
        nl = 359
    i = float(nl) / 4.
    tsky_haslam = tskylist[180*int(i) + int(j)]
    # scale temperature and add other sources of heat (50K) before returning 
    return tsky_haslam * (350/408.0)**(-2.6)+50
