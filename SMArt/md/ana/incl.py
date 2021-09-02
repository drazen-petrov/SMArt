from SMArt.incl import math, np, pd, combinations, Defaults

def normalize_RGB(*RGB, fac = 1, col_type=tuple):
    """
    normalize RGB list to 0-1
    """
    RGB = fac * np.array(RGB).mean(axis=0) / 256
    return col_type(RGB)

def fix_float(x, tolerance = 5):
    """fixes float instances using round"""
    return round(float(x), tolerance)

def get_lset(lset):
    ls = []
    for l in lset:
        ls.append(Real.fix_float(l))
    return sorted(ls)

class Real(Defaults):
    #__tolerance = 5
    #__cutoff = 10**(-5)

    def __init__(self, x):
        self.value = self.fix_float(x)

    @classmethod
    def fix_float(cls, x):
        return fix_float(x, tolerance=cls.__tolerance)

    def fix(self):
        self.value = self.fix_float(self.value)

    @classmethod
    def set_tolerance(cls, tolerance):
        assert isinstance(tolerance, int)
        cls.__tolerance = tolerance
        cls.__cutoff = 10**(-1*tolerance)

Real._add_class_defaults({'__tolerance':5, '__cutoff':10**(-5)}, flag_set=True)

kb = 0.00831451
kcal = 4.184

def _not_start_with_comm(l, comments):
    for comm in comments:
        if l.startswith(comm):
            return False
    return True

def _remove_comm(l, comments):
    for i in comments:
        l = l[:l.find(i)]
    return l

def _get_lines(f_path, comments, skip, stride):
    f = open(f_path)
    c = 0
    if skip:
        for l in f:
            if _not_start_with_comm(l, comments):
                c+=1
            if c==skip:
                break
    if stride:
        c = 0
        for l in f:
            if _not_start_with_comm(l, comments):
                if c % stride==0:
                    yield l
                c+=1
    else:
        for l in f:
            if _not_start_with_comm(l, comments):
                yield l


def read_data(f_path, skip_stride=None, **kwargs):
    if skip_stride:
        comments = kwargs.get('comments', [])
        f_path = _get_lines(f_path, comments, skip_stride[0], skip_stride[1])
    data = np.loadtxt(f_path, comments = comments, **kwargs)
    return data

# integration methods
# dgdl=[[lp1,dgdl1],[lp2,dgdl2],...] - have to be sorted

def sim(dgdl):
    """
    numerical integration using Simpson's rule
    deprecated - use scipy.integrate.simpson instead
    """
    if len(dgdl) % 2 == 0:
        print("even num of lambda points!!!")
        return
    dgdl_half = dgdl[::2]
    return (4 * trap(dgdl) - trap(dgdl_half)) / 3.

def sim_err(dgdl):
    if len(dgdl) % 2 == 0:
        print("even num of lambda points!!!")
        return
    s = (dgdl[0][1] * (4 * dgdl[1][0] - 3 * dgdl[0][0] - dgdl[2][0]))**2
    s += (dgdl[-1][1] * (3 * dgdl[-1][0] - 4 * dgdl[-2][0] + dgdl[-3][0]))**2
    for i in range(len(dgdl) - 2):
        if i % 2:
            s += (dgdl[i + 1][1] * (4 * dgdl[i + 2][0] - 4 * dgdl[i][0] - dgdl[i + 3][0] + dgdl[i - 1][0]))**2
        else:
            s += (4 * dgdl[i + 1][1] * (dgdl[i + 2][0] - dgdl[i][0]))**2
    return math.sqrt(s) / 6.


def _trap(dgdl):
    s = 0
    for i in range(len(dgdl) - 1):
        s += (dgdl[i + 1][0] - dgdl[i][0]) * (dgdl[i + 1][1] + dgdl[i][1]) / 2
    return s


def trap_err(dgdl):
    s = (dgdl[0][1] * (dgdl[1][0] - dgdl[0][0]))**2 + (dgdl[-1][1] * (dgdl[-1][0] - dgdl[-2][0]))**2
    for i in range(len(dgdl) - 2):
        s += (dgdl[i + 1][1] * (dgdl[i + 2][0] - dgdl[i][0]))**2
    return math.sqrt(s) / 2.


def trap(dgdl, interval=False):
    """
    numerical integration using trapezoidal rule
    deprecated - use scipy.integrate.trapezoid instead
    """
    if interval:
        dg = [0]
        for i in range(len(dgdl) - 1):
            dg.append(dg[-1] + _trap(dgdl[i:i + 2]))
        return dg
    else:
        return _trap(dgdl)

def CI(av1, av2, std1, std2):
    """
    Crooks intersection
    """
    if std1 == std2:
        return (av1 - av2) / 2
    if std1 == 0 or std2 == 0:
        return (av1 - av2) / 2
    temp1 = av1 / (std1**2) + av2 / (std2**2)
    temp2 = (av1 + av2)**2 / (std1**2) / (std2**2) + 2 * (1 / (std1**2) - 1 / (std2**2)) * math.log(std2 / std1)
    temp2 = math.sqrt(temp2)
    temp3 = 1 / (std1**2) - 1 / (std2**2)
    dg1 = (temp1 + temp2) / temp3
    dg2 = (temp1 - temp2) / temp3
    if av1 > (-av2):
        if dg1 < av1 and dg1 > (-av2):
            return dg1
        elif dg2 < av1 and dg2 > (-av2):
            return dg2
    else:
        if dg1 > av1 and dg1 < (-av2):
            return dg1
        elif dg2 > av1 and dg2 < (-av2):
            return dg2
    return (av1 - av2) / 2

def Jarz(dGs, T=300):
    """
    Jarzynski
    """
    beta = kb * T
    s = 0
    for i in dGs:
        s += math.exp(-i / beta)
    return -1 * beta * math.log(s / len(dGs))

def get_lifetime_trans(states_tser, N_states, state_offset=1, time_tser=None, include_first_state=True, include_last_state=True):
    """
    count number of transitions and lifetimes of each of the states
    """
    state_trans = {st:np.zeros(N_states, dtype=int) for st in range(N_states)}
    life_times = {st:[] for st in range(N_states)}
    first_state = states_tser[0]
    first_transition = None
    current_state = states_tser[0] - state_offset
    start_time = 0
    if time_tser is None:
        time_tser = range(len(states_tser))
    assert min(np.diff(time_tser)) > 0
    for t, st in zip(time_tser, states_tser):
        st -= state_offset
        if st!=current_state:
            if include_first_state:
                life_times[current_state].append(t - start_time)
                state_trans[current_state][st] += 1
            else:
                include_first_state = True
            start_time = t
            current_state = st
    if include_last_state and t!=start_time:
        life_times[current_state].append(t - start_time)
    state_trans = pd.DataFrame(state_trans)
    return life_times, state_trans

def _RMSD2(coord1, coord2):
    d2 = coord1 - coord2
    rmsd2 = (d2 * d2).sum(axis=1)
    return np.nanmean(rmsd2)

def RMSD(coord1, coord2):
    rmsd2 = _RMSD2(coord1, coord2)
    return np.sqrt(rmsd2)

def _RMSF2_coord_space(c_space):
    d2 = c_space - np.nanmean(c_space, axis = 0)
    d2 *= d2
    rmsf2 = np.nansum(d2, axis=(2,)).mean(axis = 0).mean()
    return rmsf2

def __make_c_space(coord_space, flag_make_c_space=True):
    if flag_make_c_space:
        return np.array(coord_space)
    else:
        return coord_space[0]

def _RMSD(*coord_space, flag_make_c_space=True):
    """
    simple multi conf RMSD - based on RMSF2 formula
    """
    c_space = __make_c_space(coord_space, flag_make_c_space=flag_make_c_space)
    rmsf2 = _RMSF2_coord_space(c_space)
    rmsd2 = 2 * c_space.shape[0] * rmsf2 / (c_space.shape[0] - 1)
    return np.sqrt(rmsd2)
    """
    d2 = c_space - c_space.mean(axis = 0)
    d2 *= d2
    rmsf2 = d2.sum(axis=(2,)).mean(axis = 0).mean()
    rmsd2 = 0
    for coord1, coord2 in combinations(c_space, 2):
        rmsd2 += _RMSD2(coord1, coord2)
    rmsd2 = (2 * rmsd2) / (c_space.shape[0] * (c_space.shape[0] - 1))
    np.testing.assert_almost_equal(rmsd2 * , 2 * c_space.shape[0] * rmsf2)
    return rmsd2, rmsf2
    """

def _RMSD_pairwise(*coord_space, flag_make_c_space=True):
    """
    mean square pairwise RMSD2 on multiple configurations -> normal RMSD for 2 structures
    """
    c_space = __make_c_space(coord_space, flag_make_c_space=flag_make_c_space)
    c, rmsd2 = 0, 0
    for coord1, coord2 in combinations(c_space, 2):
        temp_rmsd2 = _RMSD2(coord1, coord2)
        if not np.isnan(temp_rmsd2):
            rmsd2 += temp_rmsd2
            c += 1
    if c==0:
        return np.nan
    rmsd2 /= c
    return np.sqrt(rmsd2)

