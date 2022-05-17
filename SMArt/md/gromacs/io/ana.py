from SMArt.incl import np, pd, re
from SMArt.md.ana.incl import _not_start_with_comm, _get_lines, Real

def get_xvg_labels(f_path, comments=('#', '@')):
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

def read_xvg_data(f_path, comments=XVG_COMMENTS, skip_stride=None, **kwargs):
    if skip_stride:
        f_path = _get_lines(f_path, comments, skip_stride[0], skip_stride[1]) # generator that skips lines
    data = np.loadtxt(f_path, comments = comments, **kwargs)
    return data

def read_bar_data(f_path, dl_max=0.3, comments=XVG_COMMENTS, skip_stride=None, 
                  flag_parse_E=True, **kwargs):
    labels = get_xvg_labels(f_path, comments)
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
    if E_pos:
        lp2use_pos.append(E_pos)
    data = read_xvg_data(f_path, comments=comments, skip_stride=skip_stride, usecols=lp2use_pos, **kwargs)
    if E_pos:
        data = data[:, :-1] + data[:, -1][:, np.newaxis]
    return pd.DataFrame(data, columns = lp2use), sim_l
