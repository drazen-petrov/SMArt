from SMArt.incl import np, pd, re
from SMArt.md.ana.incl import _not_start_with_comm, _get_lines, Real

def _get_xvg_labels(f_path, comments=('#', '@')):
    f = open(f_path)
    labels = ['x']
    for l in f:
        temp_res = re.search('@[ ]*xaxis.*label[ ]*"(.*)"', l)
        if temp_res:
            labels[0] = temp_res.group(1)
        temp_res = re.search(r'@ s\d* legend[ ]*"(.*)"', l)
        if temp_res:
            labels.append(temp_res.group(1))
        if _not_start_with_comm(l, comments):
            break
    return labels

XVG_COMMENTS = ('#', '@')

def read_xvg_data(f_path, comments=XVG_COMMENTS, skip=None, stride=None, **kwargs):
    """
        parses data from a xvg file
    :param f_path: path to a xvg file
    :param comments: comment characters (skip these lines)
    :param skip: skipping first `skip` lines (int)
    :param stride: skipping every `stride` lines (int)
    :param kwargs:
        passed to np.loadtxt function
    :return:
        data - numpy array
    """
    if skip or stride:
        f_path = _get_lines(f_path, comments, skip, stride) # generator that skips lines
    data = np.loadtxt(f_path, comments = comments, **kwargs)
    return data

def read_bar_data(f_path, dl_max=0.3, comments=XVG_COMMENTS, skip=None, stride=None, 
                  flag_parse_E=True, **kwargs):
    """
        parse bar data from GROMACS *xvg file
    :param f_path: path to file
    :param dl_max: max delta lam to read (e.g. if sim_lp == 0.1 and dl_max=0.3, read data for lam in range [0., 0.4])
    :param comments: comment characters (skip these lines)
    :param skip: skipping first `skip` lines (int)
    :param stride: skipping every `stride` lines (int)
    :param flag_parse_E: parse total energies for each snapshot and add to the relative ones (needed for error estimates)
    :param kwargs:
        float_precision: number of digits used to write float in the output file - use if numbers are not separated with a space
        col_format: list of int (num of characters for each column) - use if numbers are not separated with a space.
            if only one number given, one can provide N_cols (if not, N_cols is determined from the header)
            if True given, col_format will be deduced from len(line)
        N_cols: number of columns
        also passed to read_xvg_data function
    :return:
        data as a `pandas` `DataFrame`
        simulated lambda (float)
    """
    labels = _get_xvg_labels(f_path, comments)
    lp2use = []
    lp2use_pos = []
    E_pos = None
    for i, l in enumerate(labels):
        if 'fep-lambda' in l:
            sim_l = float(re.search(r'fep-lambda.*=[ ]*(\d+(?:\.\d+)?)', l).groups()[0])
        if "\\xD\\f{}H \\xl\\f{}" in l:
            temp_l = float(l.split()[-1])
            if Real.fix_float(abs(temp_l - sim_l)) <= dl_max:
                lp2use.append(temp_l)
                lp2use_pos.append(i)
        if flag_parse_E and 'total energy' in l.lower():
            E_pos = i
    if E_pos is not None:
        lp2use_pos.append(E_pos)
    data = read_xvg_data(f_path, comments=comments, skip=skip, stride=stride, usecols=lp2use_pos, **kwargs)
    if E_pos is not None:
        data = data[:, :-1] + data[:, -1][:, np.newaxis]
    return pd.DataFrame(data, columns = lp2use), sim_l
