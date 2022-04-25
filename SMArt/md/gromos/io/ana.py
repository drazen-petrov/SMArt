from SMArt.incl import np, pd, re
from SMArt.md.ana.incl import _not_start_with_comm, _get_lines, Real

def _get_LPs2use(LPs, dl_max, sim_l, LP_offset):
    lp2use_pos = []
    for i in range(LP_offset, len(LPs)):
        if Real.fix_float(abs(LPs[i] - sim_l)) <= dl_max:
            lp2use_pos.append(i)
    return LPs[lp2use_pos], lp2use_pos

def _parse_header(f_path, dl_max, comments, N_comm_lines):
    f = open(f_path)
    c_lines_no_comm = 0
    #sim_l, step, t = None, None, None
    sim_l = None
    for c_lines, l in enumerate(f):
        if _not_start_with_comm(l, comments):
            if c_lines_no_comm==0:
                N_LPs = int(l)
            elif c_lines_no_comm==1:
                LPs = np.array(l.split(), dtype=float)
            c_lines_no_comm+=1
        #"""
        else:
            if c_lines==0 and 'sim time' in l:
                temp = l.split()
                step, t = int(temp[-2]), float(temp[-1])
            elif c_lines in (0, 1) and 'lam_s' in l:
                sim_l = float(l.split()[-1])
        #"""
        if c_lines_no_comm>N_comm_lines:
            line_len = len(l)
            break
        c_lines+=1
    LP_offset = 0
    if len(LPs)!=N_LPs:
        assert len(LPs) == N_LPs + 1
        sim_l = LPs[0]
        LP_offset = 1
    assert sim_l is not None
    lp2use, lp2use_pos = _get_LPs2use(LPs, dl_max, sim_l, LP_offset)
    return c_lines, lp2use, lp2use_pos, LPs, sim_l, len(LPs), line_len

def __get_col_pos(col_format, N_cols=None):
    col_pos = []
    if isinstance(col_format, int):
        for i in range(N_cols):
            col_pos.append((i*col_format, (i+1)*col_format))
    else:
        current_pos = 0
        for n_characters in col_format:
            if n_characters is None:
                end_pos = -1
            else:
                end_pos = current_pos + n_characters
            col_pos.append((current_pos, end_pos))
            current_pos = end_pos
    return col_pos

def __parse_line_data_collen(line_gen, col_format, N_cols=None, usecols=None):
    col_pos = __get_col_pos(col_format, N_cols)
    if usecols is None:usecols = range(len(col_pos))
    data = []
    for l in line_gen:
        temp_data = [float(l[col_pos[i][0]:col_pos[i][1]]) for i in usecols]
        data.append(temp_data)
    return np.array(data)

def __parse_line_float_precision(line_gen, float_precision=8, usecols=None):
    data = []
    for l in line_gen:
        temp_data = [float(val) for val in re.findall('[-,\d]\d*\.\d{'+str(float_precision)+'}?', l)]
        data.append(temp_data)
    data = np.array(data)
    if usecols:
        return data[:, tuple(usecols)]
    else:
        return data

def read_bar_dhdl(f_path, dl_max=0.3, comments=('#',), skip_stride=None, N_comm_lines=2, **kwargs):
    """
        parse bar / dhdl data from GROMOS (output of ext_ti_ana)
    :param f_path: path to file
    :param dl_max: max delta lam to read (e.g. if sim_lp == 0.1 and dl_max=0.3, read data for lam in range [0., 0.4])
    :param comments: comment characters (skip these lines)
    :param skip_stride: (int, int) - skipping lines
    :param N_comm_lines: number of commented lines in the header (default 2)
    :param kwargs:
        float_precision: number of digits used to write float in the output file - use if numbers are not separated with a space
        col_format: list of int (num of characters for each column) - use if numbers are not separated with a space.
            if only one number given, one can provide N_cols (if not, N_cols is determined from the header)
            if True given, col_format will be deduced from len(line)
        N_cols: number of columns
    :return: 
        data as a `pandas` `DataFrame`
        simulated lambda (float)
    """
    c_lines, lp2use, lp2use_pos, LPs, sim_l, N_cols, line_len = _parse_header(f_path, dl_max, comments, N_comm_lines)
    data = None
    line_gen = None
    if skip_stride:
        skip_stride = skip_stride[0] + N_comm_lines, skip_stride[1]
        line_gen = _get_lines(f_path, comments, *skip_stride)
    if kwargs.get('col_format'):
        col_format = kwargs.get('col_format')
        N_cols = kwargs.get('N_cols', N_cols)
        if col_format is True:
            col_format = line_len // N_cols
        if not line_gen:line_gen = _get_lines(f_path, comments, N_comm_lines, None)
        data = __parse_line_data_collen(line_gen, col_format, N_cols=N_cols, usecols=lp2use_pos)
    if kwargs.get('float_precision'):
        fl_prec = kwargs.get('float_precision')
        if not line_gen:line_gen = _get_lines(f_path, comments, N_comm_lines, None)
        data = __parse_line_float_precision(line_gen, float_precision=fl_prec, usecols=lp2use_pos)
    if data is None:
        if line_gen:
            data = np.loadtxt(line_gen, usecols=lp2use_pos, **kwargs)
        else:
            data = np.loadtxt(f_path, usecols=lp2use_pos, skiprows=c_lines, **kwargs)
    if len(data.shape) == 1:
        data = [data]
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