from SMArt.incl import math, np, Defaults

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
    print('\t', skip, stride)
    f = open(f_path)
    c = 0
    if skip:
        for l in f:
            if _not_start_with_comm(l, comments):
                c+=1
            if c==skip:
                break
    c = 0
    for l in f:
        if _not_start_with_comm(l, comments):
            if c % stride==0:
                yield l
            c+=1

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

