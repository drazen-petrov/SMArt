from SMArt.incl import sys, np, pd, Manager, Pool, integrate, bisect_left, bisect_right, Defaults, do_warn, itemgetter
from .incl import Real, get_lset
try:
    from pymbar import MBAR, timeseries
except ImportError:
    pass

from .incl import kb, kcal, CI, Jarz

class LPs_Map:
    def __init__(self, l0, l1, LPs, LPs_pred):
        self.l0 = Real(l0).value
        self.l1 = Real(l1).value
        self.LPs = LPs
        self.LPs_pred = LPs_pred
        self.LPs_interval = [] # LPs (simulated) in the interval
        self.LPs_pred_interval = [] # LPs_pred (predicted) in the interval
        self.pos_inLPs_LPs_interval = [] # position (in LPs) of LPs in the interval (before pos_LPs_interval_inLPs)
        self.pos_LPs_interval = [] # position (in LPs_pred) of LPs in the interval
        self.pos_LPs_pred_interval = [] # position (in LPs_pred) of LPs_pred in the interval
        self.pos_LPs = [] # position (in LPs_pred) of LPs (all)
        # full interval - for l0 / l1 not in LPs
        for i, lp in enumerate(LPs_pred):
            if self.l0 <= lp <= self.l1:
                self.LPs_pred_interval.append(lp)
                self.pos_LPs_pred_interval.append(i)
                if lp in LPs:
                    self.LPs_interval.append(lp)
                    self.pos_inLPs_LPs_interval.append(LPs.index(lp))
                    self.pos_LPs_interval.append(i)
                    self.pos_LPs.append(i)
            elif lp in LPs:
                self.pos_LPs.append(i)
        self.pos_inLPs_pred_interval__LPs_interval = np.array(self.pos_LPs_interval) - self.pos_LPs_interval[0]
        # position (in LPs_pred_interval) of LPs_interval (before pos_LPs_interval_inLPs_pred_interval)

        if l0 not in LPs or l1 not in LPs:
            self.pos_inLPs_LPs_interval_full = list(self.pos_inLPs_LPs_interval) # position (in LPs) of LPs in the interval
            if self.LPs_interval[0] != l0:
                self.pos_inLPs_LPs_interval_full.insert(0, self.pos_inLPs_LPs_interval[0] - 1)
            if self.LPs_interval[-1] != l1:
                self.pos_inLPs_LPs_interval_full.append(self.pos_inLPs_LPs_interval[-1] + 1)
            self.LPs_interval_full = [LPs[pos] for pos in self.pos_inLPs_LPs_interval_full] # LPs (simulated) in the interval

    def for_neigh_Es_nfr(self, Es, nfr, n_neigh=1):
        count = sum(nfr[:self.pos_LPs_interval[0]])
        for i in range(len(self.LPs_interval) - n_neigh):
            pos_LPs = self.pos_LPs_interval[i:i + 1 + n_neigh]
            temp_nfr = nfr[pos_LPs]
            fr_last = count + sum(temp_nfr)
            temp_Es = Es[pos_LPs, count:fr_last]
            yield temp_Es, temp_nfr
            count += nfr[self.pos_LPs_interval[i]]
   
    def for_LPs_Es_nfr_mul(self, Es, nfr_mul, dEs=None):
        count = 0
        for i,l in enumerate(self.LPs_pred):
            for temp_N_frms in nfr_mul[i]:
                if temp_N_frms:
                    temp_nfr = np.zeros(len(self.LPs_pred), dtype = int)
                    temp_nfr[i] = temp_N_frms
                    fr_last = count + temp_N_frms
                    temp_Es = Es[:, count:fr_last]
                    if dEs is None:
                        temp_dEs = None
                    else:
                        temp_dEs = dEs[:, count:fr_last]
                    count = fr_last
                    yield l, temp_Es, temp_dEs, temp_nfr

# combine, prep data
def combine_bar_dhdl(comb_data, bar_dhdl_data, sim_lp, flag_check_sim_lp=True, append_index=None):
    """
    combines data into a dictionary
    bar_dhdl_data is a pandas dataframe with columns LPs_pred and rows values
    """
    sim_lp = Real.fix_float(sim_lp)
    if flag_check_sim_lp:assert sim_lp in bar_dhdl_data.columns
    if sim_lp not in comb_data:
        comb_data[sim_lp] = []
    if append_index is None or append_index == len(comb_data[sim_lp]):
        comb_data[sim_lp].append(bar_dhdl_data)
    else:
        comb_data[sim_lp][append_index] = pd.concat([comb_data[sim_lp][append_index], bar_dhdl_data])


def combine_exTI(comb_data, exTI_data, sim_lp, flag_check_sim_lp=True):
    sim_lp = Real.fix_float(sim_lp)
    if flag_check_sim_lp:assert sim_lp in exTI_data.columns
    if sim_lp not in comb_data:
        comb_data[sim_lp] = exTI_data
    else:
        comb_data[sim_lp] = pd.concat([comb_data[sim_lp], exTI_data])

def get_LPs_pos(LPs_pred, LPs_list):
    LPs_pos = []
    for i,l in enumerate(LPs_list):
        if l in LPs_pred:
            LPs_pos.append(i)
    assert len(LPs_pos) == len(LPs_pred)
    return LPs_pos

def prep_mbar_input(data_bar_sys, data_dhdl_sys=None, LPs=None, LPs_pred=None, skip=1, offset=0,
                    data_frac=1., flag_bw=False, **kwargs):
    if LPs is None:
        LPs = get_lset(data_bar_sys.keys())
    LPs_pred_set = None
    for l in LPs:
        for temp_data in data_bar_sys[l]:
            if LPs_pred_set is None:
                LPs_pred_set = set(temp_data.columns)
            else:
                LPs_pred_set &= set(temp_data.columns)
        if data_dhdl_sys:
            for temp_data in data_dhdl_sys[l]:
                LPs_pred_set &= set(temp_data.columns)
    if LPs_pred is None:
        LPs_pred = get_lset(LPs_pred_set)
    else:
        assert len(LPs_pred_set - set(LPs_pred)) == 0
    for l_sim in LPs:
        assert l_sim in LPs_pred
    NFRs = {}
    Es = None
    dEs = None
    for l_sim in LPs_pred:
        if l_sim in LPs:
            NFRs[l_sim] = []
            for i in range(len(data_bar_sys[l_sim])):
                temp_data_bar = data_bar_sys[l_sim][i]
                LPs_pos = get_LPs_pos(LPs_pred, temp_data_bar.columns)
                tot_N_fr = temp_data_bar.shape[0]
                N_fr = int(tot_N_fr * data_frac)
                if flag_bw:
                    st = tot_N_fr - N_fr - offset
                    temp_end = tot_N_fr - offset
                else:
                    st = offset
                    temp_end = N_fr + offset
                temp_data = temp_data_bar.values[st:temp_end:skip, LPs_pos]
                if Es is None:
                    Es = temp_data
                else:
                    Es = np.append(Es, temp_data, axis = 0)
                NFRs[l_sim].append(temp_data.shape[0])
                if data_dhdl_sys:
                    temp_data_dhdl = data_dhdl_sys[l_sim][i]
                    assert set(temp_data_dhdl.columns[LPs_pos]) == set(temp_data_bar.columns[LPs_pos])
                    temp_data = temp_data_dhdl.values[st:temp_end:skip, LPs_pos]
                    if dEs is None:
                        dEs = temp_data
                    else:
                        dEs = np.append(dEs, temp_data, axis=0)
                    assert NFRs[l_sim][-1] == temp_data.shape[0]
    Es = Es.T
    if not dEs is None:
        dEs = dEs.T
    return Es, NFRs, dEs, LPs, LPs_pred

def get_nfr_from_NFRs(NFRs, LPs, LPs_pred):
    nfr = []
    for l in LPs_pred:
        if l in LPs:
            nfr.append(sum(NFRs[l]))
        else:
            nfr.append(0)
    return np.array(nfr)

def get_nfr_mul_from_NFRs(NFRs, LPs, LPs_pred):
    nfr_mul = []
    for l in LPs_pred:
        if l in LPs:
            nfr_mul.append(np.array(NFRs[l]))
        else:
            nfr_mul.append([0])
    return nfr_mul

def get_nfr_from_nfr_mul(nfr_mul):
    return np.array([sum(nfr_list) for nfr_list in nfr_mul])

def get_nfr_LPs(nfr):
    nfr_LPs = np.array(nfr)
    return nfr_LPs[nfr_LPs>0]

def _get_LPs_pred_exTI(exTI_data, LPs_pred = None):
    LPs_pred_set = None
    LPs = get_lset(exTI_data)
    for sim_l in exTI_data:
        if LPs_pred_set is None:
            LPs_pred_set = set(exTI_data[sim_l].columns)
        else:
            LPs_pred_set &= set(exTI_data[sim_l].columns)
    if LPs_pred is None:
        LPs_pred = get_lset(LPs_pred_set)
    else:
        assert len(LPs_pred_set - set(LPs_pred)) == 0
    return LPs, LPs_pred

def _get_exTI_data_avg(exTI_data):
    exTI_data_avg = {}
    for sim_l in exTI_data:
        exTI_data_avg[sim_l] = exTI_data[sim_l].mean()
    return exTI_data_avg

# get non-correlated data points (skip every nth step)
def si_data_bar_dhdl(data_bar_dhdl):
    si_l = {}
    for l in data_bar_dhdl:
        si_l[l] = []
        for i in range(len(data_bar_dhdl[l])):
            temp_si = timeseries.statisticalInefficiency(data_bar_dhdl[l][i].loc[:,l])
            si_l[l].append(temp_si)
    return si_l

def get_skips(si_l, LPs = None, skip = 1):
    if LPs is None:
        LPs = sorted(si_l.keys())
    skips = []
    for l in LPs:
        for temp_si in si_l[l]:
            skips.append(int(temp_si * skip) + 1)
    return skips

def si_skips_data_dEs(dEs, nfr_mul, skip = 1):
    c = 0
    skips = []
    for i in range(len(nfr_mul)):
        n_frms_list = nfr_mul[i]
        for n_frms in n_frms_list:
            if n_frms:
                c_end = c+n_frms
                temp_si = timeseries.statisticalInefficiency(dEs[i][c:c_end])
                skips.append(int(temp_si * skip) + 1)
                c = c_end
    return skips

# get reduced set of data (interval or with data frac)
def red_d_Es_nfr_mul(Es, dEs, nfr_mul, data_frac=1., offset=0, flag_bw=False, skips=1,
                    flag_rnd_offset=False, flag_bs=False, seed=None, flag_reseed=False, **kwargs):
    if (flag_rnd_offset or flag_bs) and flag_reseed:
        np.random.seed(seed)
    if isinstance(data_frac, float) or isinstance(data_frac, int) == 1:
        n_data = sum([len(i) for i in nfr_mul]) # number of simulations
        data_frac = [data_frac] * n_data
    c = 0
    c_sim = 0
    new_sample = []
    new_nfr_mul = []
    for nfr_list in nfr_mul:
        new_nfr_mul.append([])
        for temp_N_fr in nfr_list:
            if temp_N_fr:
                if isinstance(skips, int):
                    skip = skips
                else:
                    skip = skips[c_sim]
                if flag_rnd_offset:
                    offset = np.random.randint(0, skip)
                if flag_bs:
                    temp_range = list(range(c + (offset % skip), c+temp_N_fr, skip))
                    N_fr = int(len(temp_range) * data_frac[c_sim])
                    temp_range = np.random.choice(temp_range, N_fr)
                else:
                    N_fr = int(temp_N_fr * data_frac[c_sim])
                    if flag_bw:
                        end = c + temp_N_fr
                        st = end - N_fr
                    else:
                        st = c
                        end = c + N_fr
                    st += (offset % skip)
                    temp_range = list(range(st, end, skip))
                c_sim += 1
                new_sample.extend(temp_range)
                new_nfr_mul[-1].append(len(temp_range))
                c += temp_N_fr
        if not new_nfr_mul[-1]:
            new_nfr_mul[-1].append(0)
    new_Es = Es[:, new_sample]
    try:
        new_dEs = dEs[:, new_sample]
    except:
        new_dEs = None
    return new_Es, new_dEs, get_nfr_from_nfr_mul(new_nfr_mul), new_nfr_mul

# calculate dG and dG_err
def _get_mbar(Es, nfr, T):
    kT = kb * T
    return MBAR(Es/kT, nfr), kT

def calc_dg_mbar(Es, nfr, T = 300):
    mbar, kT = _get_mbar(Es, nfr, T)
    Deltaf_ij, dDeltaf_ij = mbar.getFreeEnergyDifferences()
    return Deltaf_ij*kT, dDeltaf_ij*kT, mbar

def calc_dg_bar(LPs_map, Es, nfr, T = 300):
    dg = [0]
    ddg2 = []
    OIs = []
    for temp_Es, temp_nfr in LPs_map.for_neigh_Es_nfr(Es, nfr):
        Deltaf_ij, dDeltaf_ij, mbar = calc_dg_mbar(temp_Es, temp_nfr, T=T)
        OIs.append(calc_OI_BAR(temp_Es, temp_nfr)[0])
        dg.append(dg[-1] + Deltaf_ij[0][1])
        ddg2.append(dDeltaf_ij[0][1]**2)
    return np.array(dg), np.array(ddg2), OIs

def _calc_dg_bar_LPs_pred(Es, nfr, l0, l1, LPs, LPs_pred, T = 300):
    LPs_map = LPs_Map(l0, l1, LPs, LPs_pred)
    return calc_dg_bar(LPs_map, Es, nfr, T = T)

def _calc_dg_TI(dhdl, LPs, fnc2integrate=integrate.simps):
    """calculates dG from TI"""
    dg = [0]
    for i in range(2, len(LPs) + 1):
        dg.append(fnc2integrate(dhdl[:i], LPs[:i]))
    return np.array(dg)

def calc_dg_exTI_lin(exTI_data, LPs_pred=None, flag_get_exTI_data_avg=True, fnc2integrate=integrate.simps):
    """
    calculates dG from exTI data using linear interpolation
    :param exTI_data:
    :param LPs_pred:
    :param flag_get_exTI_data_avg:
    :param fnc2call_left_right: (np.mean by default); it could also be min, max, lambda x:x[0] (return left)
    :return:
        exTI_err excluding / including simulated LPs
    """
    LPs, LPs_pred = _get_LPs_pred_exTI(exTI_data, LPs_pred)
    if flag_get_exTI_data_avg:
        exTI_data_avg =_get_exTI_data_avg(exTI_data)
    else:
        exTI_data_avg = exTI_data
    dgdl = []
    for lp in LPs_pred:
        if lp not in exTI_data:
            pos = bisect_left(LPs, lp)
            d1 = lp - LPs[pos - 1]
            d2 = LPs[pos] - lp
            v1 = exTI_data_avg[LPs[pos - 1]][lp]
            v2 = exTI_data_avg[LPs[pos]][lp]
#            dgdl_delta[l] = v2 - v1
            D_sum = d1 + d2
            new_dgdl = (v1 * d2 + v2 * d1) / D_sum
            dgdl.append(new_dgdl)
        else:
            dgdl.append(exTI_data_avg[lp][lp])
    return _calc_dg_TI(dgdl, LPs_pred, fnc2integrate = fnc2integrate)
    #return fnc2integrate(dgdl, LPs_pred)

def calc_exTI_err(exTI_data, LPs_pred = None, flag_get_exTI_data_avg=True, fnc2call_left_right = np.mean):
    """
    calculates the difference between predictions from neighbouring points
    :param exTI_data:
    :param LPs_pred:
    :param flag_get_exTI_data_avg:
    :param fnc2call_left_right: (np.mean by default); it could also be min, max, lambda x:x[0] (return left)
    :return:
        exTI_err excluding / including simulated LPs
    """
    LPs, LPs_pred = _get_LPs_pred_exTI(exTI_data, LPs_pred)
    if flag_get_exTI_data_avg:
        exTI_data_avg =_get_exTI_data_avg(exTI_data)
    else:
        exTI_data_avg = exTI_data
    dhdl_err, dhdl_err2 = [], []
    for lp in LPs_pred:
        pos = bisect_left(LPs, lp)
        if lp not in exTI_data:
            v1 = exTI_data_avg[LPs[pos - 1]][lp]
            v2 = exTI_data_avg[LPs[pos]][lp]
            err = abs(v1 - v2)
            dhdl_err.append(err)
            dhdl_err2.append(err)
        else:
            dhdl_err.append(0)
            temp_err2 = []
            v = exTI_data_avg[LPs[pos]][lp]
            if pos!=0:
                temp_err2.append(abs(v - exTI_data_avg[LPs[pos - 1]][lp]))
            if pos + 1 != len(LPs):
                temp_err2.append(abs(v - exTI_data_avg[LPs[pos + 1]][lp]))
            dhdl_err2.append(fnc2call_left_right(temp_err2))
    return np.array(dhdl_err), np.array(dhdl_err2)

def _calc_exTI_pred_mbar(mbar, dEs, LPs_pred, flag_calc_w = True):
    if flag_calc_w:
        w = mbar.getWeights().T
    else:
        w = mbar
    dhdl_pred_mbar = []
    for i in range(len(LPs_pred)):
        dhdl_pred_mbar.append(np.dot(w[i], dEs[i]))
    return np.array(dhdl_pred_mbar)

def _calc_singleLP_exTI_pred_mbar(Es, dEs, LPs_pred, nfr, T = 300):
    assert Es.shape == dEs.shape
    kT = kb * T
    mbar = MBAR(Es / kT, nfr)
    w = mbar.getWeights()
    wt = w.T
    dgdl_int = []
    for i in range(len(LPs_pred)):
        dgdl_int.append(np.dot(wt[i], dEs[i]))
    return np.array(dgdl_int)

def calc_dg_exTI_mbar(mbar, dEs, LPs_pred, flag_calc_w = True, fnc2integrate=integrate.simps):
    # calculates dG from exTI in combination with MBAR
    dhdl_pred_mbar = _calc_exTI_pred_mbar(mbar, dEs, LPs_pred, flag_calc_w=flag_calc_w)
    return _calc_dg_TI(dhdl_pred_mbar, LPs_pred, fnc2integrate = fnc2integrate), dhdl_pred_mbar
    dg = [0]
    for i in range(2, len(LPs_pred) + 1):
        dg.append(fnc2integrate(dhdl_pred_mbar[:i], LPs_pred[:i]))
    return np.array(dg), np.array(dhdl_pred_mbar)

def _get_dhdl_from_exTI_data(exTI_data, LPs=None):
    dhdl = []
    if LPs is None:
        LPs = []
    for l in sorted(exTI_data):
        if l not in LPs:
            LPs.append(l)
            dhdl.append(exTI_data[l].loc[:,l].mean())
    return LPs, dhdl

def calc_dg_TI(exTI_data):
    LPs, dhdl = _get_dhdl_from_exTI_data(exTI_data)
    return _calc_dg_TI(dhdl, LPs)

def calc_OI_BAR(Es, nfr, Nbin=100):
    # overlap integral
    dUfw_bw = [[], []]
    st_fr = 0
    for i in range(2):
        sing = {0: 1, 1: -1}[i]
        ii = (i + 1) % 2
        for j in range(nfr[i]):
            dUfw_bw[i].append((Es[ii][st_fr + j] - Es[i][st_fr + j]) * sing)
        st_fr += nfr[i]
    hist_dUfw_bw = []
    rhist = (min(min(dUfw_bw[0]), min(dUfw_bw[1])), max(max(dUfw_bw[0]), max(dUfw_bw[1])))
    for i in range(2):
        temp = np.histogram(dUfw_bw[i], Nbin, rhist)
        hist_dUfw_bw.append(temp[0] / float(nfr[i]))
    OI = 0
    for i in range(Nbin):
        temp = hist_dUfw_bw[0][i] * hist_dUfw_bw[1][i]
        if temp:
            OI += 2 * temp / (hist_dUfw_bw[0][i] + hist_dUfw_bw[1][i])
    return OI, hist_dUfw_bw

# calculate segment properties
# define segment (l0, l1, LPs, LPS_pred)
# define props to calc
#   dg_ mbar, bar, exTI_mbar, exTI_lin, TI
#       same dg from above from fw / bw on data_frac
#       BS on the same dg from abobe
#   ddg from mbar and bar
#   OI
#   exTI_err

_dg_methods = ('mbar', 'bar', 'exTI_mbar', 'exTI_lin', 'TI')
_full_dg_err = ('mbar_err', 'bar_err', 'OI', 'exTI_err')
# _data_frac_dg_err -> data_frac + fw/bw + method (mbar)
# _BS_dg_err -> bs_steps + dg_method

class dG_err_tols(Defaults):
    @staticmethod
    def _get_empty_dg_err():
        return {'full':{}, 'data_frac':{}, 'BS':{}}

    @classmethod
    def get_default_dg_err_tols(cls):
        dg_err_tols = cls._get_empty_dg_err()
        for dg_err_meth in cls.__dg_err_estimators['full']:
            if dg_err_meth == 'OI':
                tol = cls.__OI_tol
            else:
                tol = cls.__dg_tol
            dg_err_tols['full'][dg_err_meth] = tol
        for data_frac in cls.__dg_err_estimators['data_frac']:
            dg_err_tols['data_frac'][data_frac] = {}
            for fwbw in cls.__dg_err_estimators['data_frac'][data_frac]:
                dg_err_tols['data_frac'][data_frac][fwbw] = {}
                for dg_meth in cls.__dg_err_estimators['data_frac'][data_frac][fwbw]:
                    dg_err_tols['data_frac'][data_frac][fwbw][dg_meth] = cls.__dg_tol
        temp_BS_err_tols = {}
        if cls.__BS_steps:
            for dg_meth in cls.__dg_err_estimators['BS']:
                temp_BS_err_tols[dg_meth] = cls.__dg_tol
            dg_err_tols['BS'] = {cls.__BS_steps:temp_BS_err_tols}
        return dg_err_tols
    
    @staticmethod
    def _convert_dG_err2tols(dg_err_tols):
        tols = {}
        for dg_err_meth, tol in dg_err_tols['full'].items():
            tols[('full', 'segment', dg_err_meth)] = tol
        for data_frac in dg_err_tols['data_frac']:
            for fwbw in dg_err_tols['data_frac'][data_frac]:
                for dg_meth, tol in dg_err_tols['data_frac'][data_frac][fwbw].items():
                    tols[('data_frac', data_frac, fwbw, 'segment', dg_meth)] = tol
        for BS_steps in dg_err_tols['BS']:
            for dg_meth, tol in dg_err_tols['BS'][BS_steps].items():
                tols[('BS', BS_steps, 'segment', dg_meth)] = tol
        return tols

    @staticmethod
    def _convert_dG_err2keys(dg_err):
        dG_err_keys = []
        if 'full' in dg_err:
            for dg_err_meth in dg_err['full']:
                dG_err_keys.append(('full', 'segment', dg_err_meth))
        if 'BS' in dg_err:
            for BS_steps in dg_err['BS']:
                for dg_meth in dg_err['BS'][BS_steps]:
                    dG_err_keys.append(('BS', BS_steps, 'segment', dg_meth))
        return dG_err_keys

    @staticmethod
    def _get_dg_methods_data_frac(dg_err):
        dg_methods = set()
        for temp_data in dg_err['data_frac'].values():
            for dg_meth in temp_data:
                dg_methods.add(dg_meth)
        return dg_methods

_dG_err_tols_defs = {}
_dG_err_tols_defs['BS_steps'] = 100
_dG_err_tols_defs['dg_tol'] = 'kT'
_dG_err_tols_defs['OI_tol'] = 0.2

_dG_err_tols_defs['dg_err_estimators'] = {'full':('mbar_err', 'OI'), 'BS':('mbar',),
                                'data_frac':{0.5: {'fw': ('mbar',), 'bw': ('mbar',)}}}

dG_err_tols._add_class_defaults(_dG_err_tols_defs, flag_set=True)

_empty_dg_err = dG_err_tols._get_empty_dg_err()

def _generate_exTI_data(LPs_map, Es, dEs, nfr_mul, T = 300):
    exTI_data = {}
    for sim_l, temp_Es, temp_dEs, temp_nfr in LPs_map.for_LPs_Es_nfr_mul(Es, nfr_mul, dEs):
        dhdl = _calc_singleLP_exTI_pred_mbar(temp_Es, temp_dEs, LPs_map.LPs_pred, temp_nfr, T=T)
        combine_exTI(exTI_data, pd.DataFrame(dhdl.reshape((1, dhdl.shape[0])), columns=LPs_map.LPs_pred), sim_l)
    return exTI_data


def _get_dg_seg(dg, LPs_map):
    if len(dg.shape)==2:
        return dg[LPs_map.pos_LPs_pred_interval[0]][LPs_map.pos_LPs_pred_interval[-1]]
    elif dg.shape[0] == len(LPs_map.LPs_pred):
        return dg[LPs_map.pos_LPs_pred_interval[-1]] - dg[LPs_map.pos_LPs_pred_interval[0]]
    elif dg.shape[0] == len(LPs_map.LPs):
        return dg[LPs_map.pos_inLPs_LPs_interval[-1]] - dg[LPs_map.pos_inLPs_LPs_interval[0]]
    else:
        temp_text = str(dg.shape) + ', ' + str(len(LPs_map.LPs)) + ', ' + str(len(LPs_map.LPs_pred))
        raise Exception('dG shape not fitting to LPs_map: ' + temp_text)

def _get_bar_err_seg(ddg_bar2, LPs_map):
    temp_seg_pos = LPs_map.pos_inLPs_LPs_interval[0], LPs_map.pos_inLPs_LPs_interval[-1]
    return np.sqrt(ddg_bar2[temp_seg_pos[0]: temp_seg_pos[1]].sum())

def _get_OI_seg(OI, LPs_map):
    pos_inLPs_LPs_interval = getattr(LPs_map, 'pos_inLPs_LPs_interval_full', LPs_map.pos_inLPs_LPs_interval)
    temp_seg_pos = pos_inLPs_LPs_interval[0], pos_inLPs_LPs_interval[-1]
    return OI[temp_seg_pos[0]: temp_seg_pos[1]]

def _get_exTI_err_seg(exTI_err, LPs_map):
    temp_seg_pos = LPs_map.pos_LPs_pred_interval[0], LPs_map.pos_LPs_pred_interval[-1] + 1
    dhdl_err = exTI_err[temp_seg_pos[0]:temp_seg_pos[1]]
    return _calc_dg_TI(dhdl_err, LPs_map.LPs_pred_interval)[-1]

def _get_BS_err(bs_data, segs2calc_LPs_maps, temp_bs_dg):
    BS_err_lists = {}
    for seg2calc in segs2calc_LPs_maps:
        BS_err_lists[seg2calc] = {}
        for temp_meth in temp_bs_dg:
            BS_err_lists[seg2calc][temp_meth] = []
    for temp_bs_data_row in bs_data:
        temp_bs_dg_full = _get_dG_from_dG_full(temp_bs_data_row, segs2calc_LPs_maps)
        for seg2calc in segs2calc_LPs_maps:
            for temp_meth in temp_bs_dg:
                BS_err_lists[seg2calc][temp_meth].append(temp_bs_dg_full[seg2calc][temp_meth])
    BS_err = {}
    for seg2calc in segs2calc_LPs_maps:
        BS_err[seg2calc] = {}
        for temp_meth in temp_bs_dg:
            BS_err[seg2calc][temp_meth] = np.std(BS_err_lists[seg2calc][temp_meth])
    return BS_err

def _calc_seg_dG(dg_methods, LPs_map, Es, dEs, nfr_mul, T, exTI_data=None, TI_data=None):
    nfr = get_nfr_from_nfr_mul(nfr_mul)
    mbar, OI, dg_mbar, dg_bar, dg_exTI_mbar, dg_exTI_lin, dg_TI = (None, ) * 7
    dG = {}
    for dg_meth in dg_methods:
        if dg_meth == 'mbar':
            dg_mbar, ddg_mbar, mbar = calc_dg_mbar(Es, nfr, T=T)
            dG[dg_meth] = dg_mbar
        elif dg_meth == 'bar':
            dg_bar, ddg_bar2, OI = calc_dg_bar(LPs_map, Es, nfr, T=T)
            dG[dg_meth] = dg_bar
        elif dg_meth == 'exTI_mbar':
            if mbar is None:
                mbar, kT = _get_mbar(Es, nfr, T)
            dg_exTI_mbar, dhdl_pred_mbar = calc_dg_exTI_mbar(mbar, dEs, LPs_map.LPs_pred)
            dG[dg_meth] = dg_exTI_mbar
        elif dg_meth == 'exTI_lin':
            if exTI_data is None:
                exTI_data = _generate_exTI_data(LPs_map, Es, dEs, nfr_mul, T=T)
            dg_exTI_lin = calc_dg_exTI_lin(exTI_data)
            dG[dg_meth] = dg_exTI_lin
        elif dg_meth == 'TI':
            if TI_data:
                dg_TI = calc_dg_TI(TI_data)
            if exTI_data is None:
                exTI_data = _generate_exTI_data(LPs_map, Es, dEs, nfr_mul, T=T)
            dg_TI = calc_dg_TI(exTI_data)
            dG[dg_meth] = dg_TI
    return dG, mbar, OI, nfr

def _get_dG_from_dG_full(dG_full, segs2calc_LPs_maps):
    dG = dict((seg2calc, {}) for seg2calc in segs2calc_LPs_maps)
    for dg_meth in dG_full:
        for seg2calc, temp_LPs_map in segs2calc_LPs_maps.items():
            dG[seg2calc][dg_meth] = _get_dg_seg(dG_full[dg_meth], temp_LPs_map)
    return dG

def calc_seg_props(LPs_map, Es, dEs, nfr_mul, si_skips=None, T=300, **kwargs):
    """
    flag_calc_full_seg_props - include full segment in calc_props
    """
    segs2calc_LPs_maps = kwargs.get('segs2calc_LPs_maps', {})
    if not segs2calc_LPs_maps:
        if 'segs2calc' in kwargs:
            segs2calc = kwargs.get('segs2calc')
            for seg2calc in segs2calc: 
                segs2calc_LPs_maps[seg2calc] = LPs_Map(seg2calc[0], seg2calc[1], LPs_map.LPs, LPs_map.LPs_pred)
        else:
            segs2calc = [(LPs_map.l0, LPs_map.l1)]
            segs2calc_LPs_maps[(LPs_map.l0, LPs_map.l1)] = LPs_map
    if kwargs.get('flag_calc_full_seg_props'):
        segs2calc_LPs_maps[(LPs_map.l0, LPs_map.l1)] = LPs_map
         # this is nice if seg = (0., 1.) - which gives than the full dG and dG_err
    dg_methods = kwargs.get('dg', {'mbar'})
    dg_err_estimators = kwargs.get('dg_err', dG_err_tols.get_default_dg_err_tols())
    dg_methods |= dG_err_tols._get_dg_methods_data_frac(dg_err_estimators)
    exTI_data = kwargs.get('exTI_data')
    TI_data = kwargs.get('TI_data')
    dG_full, mbar, OI, nfr = _calc_seg_dG(dg_methods, LPs_map, Es, dEs, nfr_mul, T, exTI_data, TI_data)
    # dG calculations
    dG = _get_dG_from_dG_full(dG_full, segs2calc_LPs_maps)
    # dG_err calculations
    dG_err = {'full':{}, 'data_frac':{}, 'BS':{}}
    for seg2calc in segs2calc_LPs_maps:
        dG_err['full'][seg2calc] = {}
    exTI_err, exTI_err2 = None, None
    full_dg_err_list = dg_err_estimators.get('full',[])
    if {'mbar_err', 'bar_err'} & set(full_dg_err_list) or dg_err_estimators['BS']:
        if si_skips is None:
            si_skips = si_skips_data_dEs(Es, nfr_mul)
        elif isinstance(si_skips, dict):
            si_skips = get_skips(si_skips, LPs_map.LPs)
        si_Es, si_dEs, si_nfr, si_nfr_mul = red_d_Es_nfr_mul(Es, None, nfr_mul, skips = si_skips)
    for full_dg_err in dg_err_estimators.get('full', []):
        if full_dg_err == 'mbar_err':
            si_dg_mbar, si_ddg_mbar, si_mbar = calc_dg_mbar(si_Es, si_nfr, T=T)
            for seg2calc, temp_LPs_map in segs2calc_LPs_maps.items():
                dG_err['full'][seg2calc][full_dg_err] = _get_dg_seg(si_ddg_mbar, temp_LPs_map)
#            dG_err['full'][full_dg_err] = si_ddg_mbar
        elif full_dg_err == 'bar_err':
            si_dg_bar, si_ddg_bar2, si_OI = calc_dg_bar(LPs_map, si_Es, si_nfr, T=T)
            for seg2calc, temp_LPs_map in segs2calc_LPs_maps.items():
                dG_err['full'][seg2calc][full_dg_err] = _get_bar_err_seg(si_ddg_bar2, temp_LPs_map)
#            dG_err['full'][full_dg_err] = si_ddg_bar2
        elif full_dg_err == 'OI':
            if OI is None:
                dg_bar, ddg_bar2, OI = calc_dg_bar(LPs_map, Es, nfr, T=T)
            for seg2calc, temp_LPs_map in segs2calc_LPs_maps.items():
                dG_err['full'][seg2calc][full_dg_err] = _get_OI_seg(OI, temp_LPs_map)
#            dG_err['full'][full_dg_err] = OI
        elif full_dg_err == 'exTI_err':
            if exTI_data is None:
                exTI_data = _generate_exTI_data(LPs_map, Es, dEs, nfr_mul, T=T)
            exTI_err, exTI_err2 = calc_exTI_err(exTI_data)
            for seg2calc, temp_LPs_map in segs2calc_LPs_maps.items():
                dG_err['full'][seg2calc][full_dg_err] = _get_exTI_err_seg(exTI_err, temp_LPs_map)
                dG_err['full'][seg2calc][full_dg_err + '2'] = _get_exTI_err_seg(exTI_err2, temp_LPs_map)
#            dG_err['full'][full_dg_err] = exTI_err
#            dG_err['full'][full_dg_err + '2'] = exTI_err2
    # data fraction
    data_frac_dG = dG_err['data_frac']
    data_frac_dg_err_list = dg_err_estimators.get('data_frac', {})
    for data_frac in data_frac_dg_err_list:
        data_frac_dG[data_frac] = {}
        for fwbw in data_frac_dg_err_list[data_frac]:
            assert fwbw in ('fw', 'bw')
            if fwbw == 'fw':
                flag_bw = False
            else:
                flag_bw = True
            temp_dEs = None
            temp_dg_methods = set(data_frac_dg_err_list[data_frac][fwbw])
            if set(('exTI_mbar', 'exTI_lin', 'TI')) & temp_dg_methods:
                temp_dEs = dEs
            temp_Es, temp_dEs, temp_nfr, temp_nfr_mul = red_d_Es_nfr_mul(Es, dEs, nfr_mul, data_frac=data_frac, flag_bw=flag_bw)
            data_frac_dG_full= _calc_seg_dG(temp_dg_methods, LPs_map, temp_Es, temp_dEs, temp_nfr_mul, T)[0]
            temp_data_frac_dG = _get_dG_from_dG_full(data_frac_dG_full, segs2calc_LPs_maps)
            for temp_seg in temp_data_frac_dG:
                for temp_dg_meth in temp_data_frac_dG[temp_seg]:
                    temp_diff = abs(temp_data_frac_dG[temp_seg][temp_dg_meth] - dG[temp_seg][temp_dg_meth])
                    temp_data_frac_dG[temp_seg][temp_dg_meth] = temp_diff
            data_frac_dG[data_frac][fwbw] = temp_data_frac_dG
    # bootstrap err
    if dg_err_estimators['BS']:
        for n_bs_steps in dg_err_estimators['BS']:
            temp_bs_dg = dg_err_estimators['BS'][n_bs_steps]
            ncpu = kwargs.get('ncpu')
            bs_kwargs = {}
            if 'seed_MAX' in kwargs:
                bs_kwargs['seed_MAX'] = kwargs['seed_MAX']
            if ncpu:
                bs_data = _do_bs_ana_Ncpu(n_bs_steps, temp_bs_dg, LPs_map, Es, dEs, nfr_mul, si_skips, T, ncpu=ncpu, **bs_kwargs)
            else:
                #bs_data = _do_bs_ana_1cpu(n_bs_steps, temp_bs_dg, LPs_map, Es, dEs, nfr_mul, si_skips, T)
                bs_data = _do_bs_ana_Ncpu(n_bs_steps, temp_bs_dg, LPs_map, Es, dEs, nfr_mul, si_skips, T, ncpu=1, **bs_kwargs)
            dG_err['BS'][n_bs_steps] = _get_BS_err(bs_data, segs2calc_LPs_maps, temp_bs_dg)
    return dG, dG_err, dG_full

"""
def _do_bs_step(args):
    dg, LPs_map, Es, dEs, nfr_mul, si_skips, T, bs_data = args
    temp_Es, temp_dEs, temp_nfr, temp_nfr_mul = red_d_Es_nfr_mul(Es, dEs, nfr_mul, skips=si_skips, flag_bs=True, flag_rnd_offset=True)
    bs_data.append(calc_seg_props(LPs_map, temp_Es, temp_dEs, temp_nfr_mul, T=T, dg=dg, dg_err=_empty_dg_err)[0])

def _do_bs_ana_1cpu(bs_steps, dg, LPs_map, Es, dEs, nfr_mul, si_skips, T, **kwargs):
    raw_bs_data = []
    for i in range(bs_steps):
        _do_bs_step((dg, LPs_map, Es, dEs, nfr_mul, si_skips, T, raw_bs_data))
    return raw_bs_data
"""

def _do_bs_step_mp(args):
    ns, temp_seed, LPs_map = args
    dg_methods, Es, dEs, nfr_mul, si_skips, T = ns.dg, ns.Es, ns.dEs, ns.nfr_mul, ns.si_skips, ns.T
    temp_Es, temp_dEs, temp_nfr, temp_nfr_mul = red_d_Es_nfr_mul(Es, dEs, nfr_mul, skips=si_skips, flag_bs=True, flag_rnd_offset=True, seed=temp_seed, flag_reseed=True)
    ns.raw_bs_data.append(_calc_seg_dG(dg_methods, LPs_map, temp_Es, temp_dEs, temp_nfr_mul, T)[0])

def _do_bs_ana_Ncpu(bs_steps, dg, LPs_map, Es, dEs, nfr_mul, si_skips, T, ncpu=4, seed_MAX=2**31, **kwargs):
    p = Pool(ncpu)
    try:
        man = Manager()
        ns = man.Namespace()
        raw_bs_data = man.list()
        ns.raw_bs_data = raw_bs_data
        #ns.raw_bs_data = man.list()
        ns.dg, ns.Es, ns.dEs, ns.nfr_mul, ns.si_skips, ns.T = dg, Es, dEs, nfr_mul, si_skips, T
        mp_args = []
        seeds = set()
        for i in range(bs_steps):
            temp_seed = np.random.randint(seed_MAX)
            while temp_seed in seeds:
                temp_seed = np.random.randint(seed_MAX)
            seeds.add(temp_seed)
            mp_args.append((ns, temp_seed, LPs_map))
        p.map_async(_do_bs_step_mp, mp_args).get()
    except:
        print("Unexpected error:\n")
        print((sys.exc_info()))
    p.terminate()
    return raw_bs_data

# update schemes
def _get_merge_segments(LPs_allowed, dl_merge_list, slide_win = None, **kwargs):
    segs2check = set()
    for dl_merge in dl_merge_list:
        if slide_win:
            temp_step = slide_win
        else:
            temp_step = dl_merge
        LPs = get_lset(np.arange(0., 1.0001 - temp_step, temp_step))
        for l in LPs:
            temp_seg = (l, Real.fix_float(l + dl_merge))
            if temp_seg[0] in LPs_allowed and temp_seg[1] in LPs_allowed:
                segs2check.add(temp_seg)
    return segs2check

def _get_additional_seg2calc(LPs, data_max_dl, segs2check):
    for temp_seg in segs2check:
        pass

def _get_LPs_in_seg(LPs, temp_seg):
    pos = bisect_left(LPs, temp_seg[0])
    pos2 = bisect_right(LPs, temp_seg[1])
    return tuple(LPs[pos:pos2])

def _add_midpoints(LPs, dl_min):
    LPs_mid = get_lset(np.array(LPs)[:-1] + np.diff(LPs) * 0.5)
    LPs_mid = set(LPs_mid) & set(get_lset(np.arange(0.,1.0000001, dl_min)))
    new_LPs = list(LPs)
    new_LPs.extend(LPs_mid)
    return get_lset(new_LPs)

def get_segments2test(LPs, LPs_allowed, seg_width_slide = [(0.3, 0.1), (0.2, 0.2)], dl_merge_list = [0.2, 0.1, 0.05], **kwargs):
    """
    :param LPs: 
    :param seg_width_slide: parameters to generate segments to calculate properties
    :param dl_merge_list: additional segments to take into account in addition ones defined by the LPs
    :param flag_midpoints: add midpoints between LPs
    :param kwargs: 
        slide_win - in combination with dl_merge_list
    :return: 
    """
    max_dl_LPs = np.diff(LPs).max()
    if dl_merge_list:
        max_dl_LPs = max(max_dl_LPs, max(dl_merge_list))
    max_dl_LPs = Real.fix_float(max_dl_LPs)
    LPs_set = set(LPs_allowed)
    segs2check = _get_merge_segments(LPs_set, dl_merge_list, **kwargs)
    for i in range(len(LPs_allowed) - 1):
        segs2check.add((LPs_allowed[i], LPs_allowed[i+1]))
    segs2calc = {}
    segs2calc_dG = []
    for seg_width, seg_slide_step in seg_width_slide:
        temp_start_LPs = get_lset(np.arange(0., 1.0001 - seg_width, seg_slide_step))
        lim_val = (seg_width - seg_slide_step) / 2 # lim value added on both sides of the segment
        cp_segs2check = set(segs2check)
        temp_segscalc_dG = {}
        for l in temp_start_LPs:
            temp_seg = (l, Real.fix_float(l + seg_width))
            if l==0.:
                temp_seg_lim = [l]
            else:
                temp_seg_lim = [Real.fix_float(l + lim_val)]
            if temp_seg[1] == 1.:
                temp_seg_lim.append(1.)
            else:
                temp_seg_lim.append(Real.fix_float(temp_seg[1] - lim_val))
            temp_seg_lim = tuple(temp_seg_lim)
            segs2calc[temp_seg] = []
            for temp_seg2check in segs2check:
                if temp_seg2check[0] >= temp_seg_lim[0] and temp_seg2check[1] <= temp_seg_lim[1]:
                    segs2calc[temp_seg].append(temp_seg2check)
            segs2check -= set(segs2calc[temp_seg])
            # save the temp_seg_lim for each temp_LPs (to calculate full dG from the segments)
            temp_LPs = _get_LPs_in_seg(LPs, temp_seg)
            temp_segscalc_dG[temp_LPs] = temp_seg_lim
        segs2calc_dG.append(temp_segscalc_dG)
    assert temp_seg[1] == 1.
    assert len(segs2check) == 0
    segs2calc_LPs = {}
    for temp_seg in segs2calc:
        if len(segs2calc[temp_seg])!=0:
            temp_LPs = _get_LPs_in_seg(LPs, temp_seg)
            if temp_LPs not in segs2calc_LPs:
                segs2calc_LPs[temp_LPs] = segs2calc[temp_seg]
            else:
                segs2calc_LPs[temp_LPs].extend(segs2calc[temp_seg])
    return segs2calc_LPs, segs2calc_dG

def __get_dG_err_data_from_keys(dG_err, dG_err_keys, temp_seg):
    flag_OI = False
    temp_data = dG_err
    while dG_err_keys:
        temp_key = dG_err_keys.pop()
        if temp_key == 'segment':
            temp_key = temp_seg
        if temp_key == 'OI':
            flag_OI = True
        temp_data = temp_data[temp_key]
    return temp_data, flag_OI

def seg_check_conv_calc_score(temp_seg, dG_err, tols, fnc_OI=min):
    score = []
    flag_conv = []
    temp_dl = Real.fix_float(temp_seg[1] - temp_seg[0])
    for tol in tols:
        flag_OI = False
        tol_keys = list(reversed(tol))
        temp_data, flag_OI = __get_dG_err_data_from_keys(dG_err, tol_keys, temp_seg)
        if flag_OI:
            temp_data = fnc_OI(temp_data)
            score.append(1 - temp_data)
            if tols[tol] is not None:
                flag_conv.append(temp_data >= tols[tol])
        else:
            score.append(temp_data)
            if tols[tol] is not None:
                flag_conv.append(temp_data <= tols[tol] * temp_dl)
    return score, sum(flag_conv) == len(flag_conv)

def check_overlapping_segs(seg1, seg2):
    seg1 = get_lset(seg1)
    seg2 = get_lset(seg2)
    if seg1[0]>seg2[0]:
        if seg1[0]>=seg2[0] and seg1[0]<=seg2[1]:
            return True
    else:
        if seg2[0]>=seg1[0] and seg2[0]<=seg1[1]:
            return True
    return False

def update_overlapping_segs(seg, segs, flag_add_nonoverlapping = True):
    for i in range(len(segs)):
        if check_overlapping_segs(seg, segs[i]):
            temp_seg = segs.pop(i)
            new_seg = (min(seg[0], temp_seg[0]), max(seg[1], temp_seg[1]))
            temp_flag = update_overlapping_segs(new_seg, segs)
            if not temp_flag and not flag_add_nonoverlapping:
                segs.append(new_seg)
            return True
    if flag_add_nonoverlapping:
        segs.append(seg)
    return False

def reset_overlapping_segs(segs):
    new_segs = [tuple(segs[0])]
    for i in range(1,len(segs)):
        update_overlapping_segs(tuple(segs[i]), new_segs)
    return new_segs

def check_fullseg_in_segs(seg, segs):
    for temp_seg in segs:
        if seg[0]>=temp_seg[0] and seg[1]<=temp_seg[1]:
            return True
    return False

def update_LPs_times(data_bar_sys, data_dhdl_sys=None, T=300, **kwargs):
    """
    :param data_bar_sys: 
    :param data_dhdl_sys: 
    :param T: default 300
    :param kwargs: 
        dg_err_tols (tolerance for different error estimates - dG_err_tols.get_default_dg_err_tols())
        seg_width_slide
        dl_merge_list
        midpoints_dl_min
    :return: 
    """
    kT = kb * T
    #data_dl_max_list = kwargs.get('data_dl_max_list', [(0.3, 0.1), (0.2, 0.2)])
    #dl_merge_list = kwargs.get('dl_merge_list', [0.2, 0.1, 0.05])
    dg_err_tols = kwargs.get('dg_err_tols', dG_err_tols.get_default_dg_err_tols())
    tols = dG_err_tols._convert_dG_err2tols(dg_err_tols)
    dl_min = Real.fix_float(kwargs.get('dl_min', 0.025))
    update_LPs_fnc2call = kwargs.get('update_LPs_fnc2call', 2)
    for tol_key in tols:
        if tols[tol_key] == 'kT':
            tols[tol_key] = kT
        elif tols[tol_key] == 'half_kT':
            tols[tol_key] = kT * 0.5
        elif tols[tol_key] == 'kcal':
            tols[tol_key] = kcal
    seg_test_kwargs = dict(kwargs)
    seg_test_kwargs['dg_err'] = dg_err_tols
    LPs_sim = get_lset(data_bar_sys)
    LPs = LPs_sim
    dl_min_update = 100000 # this basically means no midpoints in the update (if they are added in LPs_allowed)
    LPs_allowed = kwargs.get('LPs_allowed')
    if LPs_allowed is None:
        if kwargs.get('add_midpoints_LPs_allowed', True):
            LPs_allowed = _add_midpoints(LPs_sim, dl_min)
        else:
            LPs_allowed = list(LPs_sim)
            dl_min_update = dl_min # add midpoints at the update
    segs2calc_LPs, segs2calc_dG = get_segments2test(LPs_sim, LPs_allowed, **seg_test_kwargs)
    si_l_dict = si_data_bar_dhdl(data_bar_sys)
    seg_score_flag = {}
    converged_segments = []
    seg_data_dG_err = {}
    for temp_LPs, segs2calc in segs2calc_LPs.items():
        Es, NFRs, dEs, LPs, LPs_pred = prep_mbar_input(data_bar_sys, data_dhdl_sys, LPs = temp_LPs, **kwargs)
        nfr_mul = get_nfr_mul_from_NFRs(NFRs, LPs, LPs_pred)
        LPs_map = LPs_Map(LPs[0], LPs[-1], LPs, LPs_pred)
        seg_dG, seg_dG_err, dG_full = calc_seg_props(LPs_map, Es, dEs, nfr_mul, segs2calc = segs2calc, si_skips = si_l_dict, **kwargs)
        seg_data_dG_err[temp_LPs] = seg_dG, seg_dG_err, dG_full
        for temp_seg in segs2calc:
            score, flag_conv = seg_check_conv_calc_score(temp_seg, seg_dG_err, tols)
            seg_score_flag[temp_seg] = score, flag_conv
            if flag_conv:
                update_overlapping_segs(temp_seg, converged_segments)
    prev_step_converged_segments = kwargs.get('converged_segments')
    if prev_step_converged_segments:
        for temp_seg in prev_step_converged_segments:
            if not check_fullseg_in_segs(temp_seg, converged_segments):
                temp_txt = 'non-converged segments in previously converged ones:\nprevious: '
                do_warn(temp_txt + str(prev_step_converged_segments) + '\ncurrent: ' + str(converged_segments))
                break
        for temp_seg in prev_step_converged_segments:
            update_overlapping_segs(temp_seg, converged_segments)
    new_LPs = update_LPs_call(update_LPs_fnc2call, LPs_allowed, seg_score_flag, converged_segments, dl_min_update)
    return seg_score_flag, converged_segments, new_LPs, seg_data_dG_err, segs2calc_dG

def update_LPs_call(fnc2call, LPs, seg_score_flag, converged_segments, dl_min):
    update_LPs_fnc_dict = {1:update_LPs_1, 2:update_LPs_2}
    return update_LPs_fnc_dict[fnc2call](LPs, seg_score_flag, converged_segments, dl_min)

def update_LPs_1(LPs, seg_score_flag, converged_segments, dl_min):
    new_LPs = []
    for i in range(len(LPs) - 1):
        temp_seg = LPs[i], LPs[i+1]
        if not check_fullseg_in_segs(temp_seg, converged_segments):
            new_LPs.append(temp_seg[0])
            temp_new_lp = Real.fix_float(np.mean(temp_seg))
            if dl_min <= Real.fix_float(temp_seg[1] - temp_new_lp):
                new_LPs.append(temp_new_lp)
            new_LPs.append(temp_seg[1])
    return new_LPs

def update_LPs_2(LPs, seg_score_flag, converged_segments, dl_min):
    segs2sim = []
    W_segs = []
    for i in range(len(LPs) - 1):
        temp_seg = LPs[i], LPs[i+1]
        if not check_fullseg_in_segs(temp_seg, converged_segments):
            W_segs.append(seg_score_flag[temp_seg][0])
            segs2sim.append(temp_seg)
    W_segs = np.array(W_segs).T
    for i in range(len(W_segs)):
        W_segs[i]/=float(sum(W_segs[i]))
    new_LPs_weights = {}
    for i in range(len(segs2sim)):
        temp_seg = segs2sim[i]
        LPs_2_sim = [temp_seg[0]]
        temp_new_lp = Real.fix_float(np.mean(temp_seg))
        if dl_min <= Real.fix_float(temp_seg[1] - temp_new_lp):
            LPs_2_sim.append(temp_new_lp)
        LPs_2_sim.append(temp_seg[1])
        temp_w = np.average(W_segs.T[i]) / len(LPs_2_sim)
        for temp_l in LPs_2_sim:
            if temp_l not in new_LPs_weights:
                new_LPs_weights[temp_l] = 0
            new_LPs_weights[temp_l] += temp_w
    return new_LPs_weights

def get_LPs_times(new_LPs_weights, LPs_times, max_iter_LPs_t, max_iter_t_LP=1., max_t_LP=5., max_total_LPs_t=100., min_t_LP=0.5, t_step=0.1):
    """
    calculates (suggestions) simulation time (in ns) for all new lamdba points
    :param new_LPs_weights: dict of new lambda points (keys) and weights for calulations of sim time (values)
    :param LPs_times: dict of all LPs and simulation times (so far)
    :param max_iter_LPs_t: max total time for the next iteration (sum of times for all new_LPs)
    :param max_iter_t_LP: max time for each of the new LPs (in this iteration)
    :param max_t_LP: max total time for each of the new LPs (this iteration + from LPs_times)
    :param max_total_LPs_t: max total time for all LPs ad for all iterations
    :param min_t_LP: min time for each new lambda point (not in LPs_times)
    :param t_step: simulation time for new LPs is added in steps
    """
    iter_LPs_times = dict((temp_l, 0) for temp_l in new_LPs_weights)
    t_step = float(t_step)
    weight_step = t_step / max_iter_LPs_t
    if min_t_LP>t_step:
        for temp_l in new_LPs_weights:
            if temp_l not in LPs_times:
                LPs_times[temp_l] = min_t_LP
                iter_LPs_times[temp_l] += min_t_LP
                new_LPs_weights[temp_l] -= weight_step * min_t_LP / t_step
                max_iter_LPs_t -= min_t_LP
    new_LPs_weights = dict(new_LPs_weights)
    while Real.fix_float(sum(LPs_times.values()) + t_step) < max_total_LPs_t and new_LPs_weights and max_iter_LPs_t >= t_step:
        temp_l = max(new_LPs_weights.items(), key=lambda x:x[1])[0]
        if temp_l not in LPs_times:
            new_value, new_value_iter = t_step, t_step
        else:
            new_value = Real.fix_float(LPs_times[temp_l] + t_step)
            new_value_iter = Real.fix_float(iter_LPs_times[temp_l] + t_step)
        if new_value > max_t_LP or new_value_iter > max_iter_t_LP:
            new_LPs_weights.pop(temp_l)
        else:
            LPs_times[temp_l] = new_value
            iter_LPs_times[temp_l] = new_value_iter
            new_LPs_weights[temp_l] -= weight_step
            max_iter_LPs_t -= t_step
            Real.fix_float(max_iter_LPs_t)
    return iter_LPs_times


def __get_max_dl(segs):
    dls = [seg[-1] - seg[0] for seg in segs]
    return max(dls)

def get_full_dG_from_segs(seg_data_dG_err, segs2calc_dG, method='mbar'):
    """
    get the full dG from the segments
    :param seg_data_dG_err: output from update_LPs_times
    :param segs2calc_dG: output from update_LPs_times
    :param method: method used to calculate the dG of the segments, e.g. 'mbar'
    :return: dG
    """
    dG = 0
    for LPs_2_calc, seg2calc_dG in segs2calc_dG.items():
        dG += seg_data_dG_err[LPs_2_calc][0][seg2calc_dG][method]
    return dG

def get_full_dG_err_from_segs(seg_data_dG_err, segs2calc_dG, err_method=dict(full=['mbar_err'])):
    """
    get the full dG from the segments
    :param seg_data_dG_err: output from update_LPs_times
    :param segs2calc_dG: output from update_LPs_times
    :param err_method: method used to calculate the error estimate of the segments, e.g. dict(full=['mbar_err'], BS={N_steps:['mbar']})
    :return: dG
    """
    dG_err = []
    dG_err_keys = dG_err_tols._convert_dG_err2keys(err_method)
    for dG_err_key in dG_err_keys:
        temp_dG_err = 0
        for LPs_2_calc, seg2calc_dG in segs2calc_dG.items():
            temp_key = list(reversed(dG_err_key))
            temp_data, flag_OI = __get_dG_err_data_from_keys(seg_data_dG_err[LPs_2_calc][1], temp_key, seg2calc_dG)
            temp_dG_err += temp_data * temp_data
        dG_err.append(np.sqrt(temp_dG_err))
    return dG_err


