from SMArt.incl import np, pd
from SMArt.md.ana.incl import _not_start_with_comm, _get_lines, Real

def _get_LPs2use(LPs, dl_max, sim_l, LP_offset):
    lp2use_pos = []
    for i in range(LP_offset, len(LPs)):
        if Real.fix_float(abs(LPs[i] - sim_l)) <= dl_max:
            lp2use_pos.append(i)
    return LPs[lp2use_pos], lp2use_pos

def _parse_header(f_path, dl_max, comments):
    f = open(f_path)
    c_lines_no_comm, c_lines = 0, 0
    sim_l, step, t = None, None, None
    for l in f:
        if _not_start_with_comm(l, comments):
            if c_lines_no_comm==0:
                N_LPs = int(l)
            elif c_lines_no_comm==1:
                LPs = np.array(l.split(), dtype=float)
            c_lines_no_comm+=1
        else:
            if c_lines==0:
                temp = l.split()
                step, t = int(temp[-2]), float(temp[-1])
            elif c_lines==1 and 'lam_s' in l:
                sim_l = float(l.split()[-1])
        if c_lines_no_comm==3:
            break
        c_lines+=1
    LP_offset = 0
    if len(LPs)!=N_LPs:
        assert len(LPs) == N_LPs + 1
        sim_l = LPs[0]
        LP_offset = 1
    lp2use, lp2use_pos = _get_LPs2use(LPs, dl_max, sim_l, LP_offset)
    return c_lines, lp2use, lp2use_pos, LPs, sim_l, step, t

def read_bar_dhdl(f_path, dl_max=0.3, comments=('#',), skip_stride=None, **kwargs):
    c_lines, lp2use, lp2use_pos, LPs, sim_l, step, t = _parse_header(f_path, dl_max, comments)
    if skip_stride:
        skip_stride = skip_stride[0] + 2, skip_stride[1]
        line_gen = _get_lines(f_path, comments, *skip_stride)
        data = np.loadtxt(line_gen, usecols=lp2use_pos, **kwargs)
    else:
        data = np.loadtxt(f_path, usecols=lp2use_pos, skiprows=c_lines, **kwargs)
    return pd.DataFrame(data, columns = lp2use), sim_l

def get_sim_l_from_fpath(f_path):
    with open(f_path) as f:
        return float(f.readline().split()[1])

def read_exTI(f_path, sim_l = None, fnc2call=get_sim_l_from_fpath):
    f = open(f_path)
    if sim_l is None:
        sim_l = fnc2call(f_path)
    temp_data = np.loadtxt(f_path)
    return sim_l, pd.DataFrame([temp_data.T[1]], columns = temp_data.T[0])