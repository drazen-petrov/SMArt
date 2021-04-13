import random
from SMArt.md.gromos.io.incl import GromosFile


def _make_job_list(f_path, params, title='job list', flag_subdir=False):
    gf = GromosFile(f_path, write=1)
    gf.write_block('TITLE', title)
    n_sims = len(params[list(params.keys())[0]])
    s = ''
    l = {}
    if 'subdir' not in params:
        params['subdir'] = ['.'] * n_sims
    if 'job_id' not in params:
        params['job_id'] = list(range(1, n_sims + 1))
    if 'run_after' not in params:
        params['run_after'] = [0]
        params['run_after'].extend(params['job_id'][:-1])
    param_keys = list(params.keys())
    for p in ('job_id', 'run_after', 'subdir'):
        param_keys.remove(p)
    param_keys.insert(0, 'job_id')
    param_keys.append('subdir')
    param_keys.append('run_after')
    for p in param_keys:
        if len(p) > 9:
            temp = '{:>' + str(len(p) + 1) + '}'
            s += temp.format(p)
            l[p] = str(len(p) + 1)
        else:
            s += '{:>10}'.format(p)
            l[p] = '10'
    s += '\n'
    for i in range(n_sims):
        for p in param_keys:
            temp = '{:>' + l[p] + '}'
            s += temp.format(params[p][i])
        s += '\n'
    gf.write_block('JOBSCRIPTS', s)
    gf.f.close()

make_job_list = _make_job_list


def _make_T_params(T_list, n_grp):
    temp_d = {'TEMPI': T_list}
    for i in range(n_grp):
        temp_d['TEMP0[' + str(i + 1) + ']'] = T_list
    return temp_d


def make_eq_job_list(f_path, use_posres=False, use_rot_trans=False, **kwargs):
    T = list(range(50, 301, 50))
    T.append(300)
    press = [0] * (len(T) - 1)
    press.append(2)
    gen_v = [1]
    gen_v.extend([0] * (len(T) - 1))
    nsteps = [25000] * (len(T) - 1)
    nsteps.append(100000)
    rnd = [random.randint(1, 1000000) for _ in range(len(T))]
    params = _make_T_params(T, 2)
    params['NTIVEL'] = gen_v
    params['NSTLIM'] = nsteps
    params['IG'] = rnd
    params['COUPLE'] = press
    if use_posres:
        temp_cpor = 25000
        params['NTPOR'] = []
        params['CPOR'] = []
        for i in range(len(T)):
            params['CPOR'].append(str(temp_cpor))
            params['NTPOR'].append('1')
            temp_cpor /= 10
    if use_rot_trans:
        params['NSCM'] = [1000]
        params['RTC'] = [0]
        for i in range(1, len(T)):
            if params['CPOR'][i] > 0:
                params['NSCM'].append(1000)
                params['RTC'].append(0)
            else:
                params['NSCM'].append(0)
                params['RTC'].append(1)
    _make_job_list(f_path, params)
