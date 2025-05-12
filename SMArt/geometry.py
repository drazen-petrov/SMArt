from SMArt.incl import np
from SMArt.incl import scipy
if scipy:
    _align_vectors = scipy.spatial.transform.Rotation.align_vectors

def rot_2D(v, theta = 90):
    theta = np.radians(theta)
    sin = np.sin(theta)
    cos = np.cos(theta)
    R = np.array([[cos, -sin], [sin, cos]])
    return R.dot(v)

def rot(v, theta=90, axis=None):
    if axis is None:
        axis = [0, 0, 1]
    axis = np.asarray(axis)
    theta = np.radians(theta) / 2.0
    axis = axis / np.linalg.norm(axis)
    a = np.cos(theta)
    b, c, d = -axis * np.sin(theta)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    m = np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                  [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                  [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
    return np.dot(m, v)

def add_ring_path_simple(path2add, coord_df, anchor_point, temp_v, new_ap, v_fact, **kwargs):
    if len(anchor_point) == 2:
        perp_v = np.diff(coord_df.loc[list(anchor_point)].values.T)[:,0]
        perp_v /= np.linalg.norm(perp_v)
    else:
        perp_v = rot_2D(temp_v)
    temp_f = np.linspace(-0.5, 0.5, len(path2add[anchor_point]))
    for i in range(len(path2add[anchor_point])):
        coord_df.loc[path2add[anchor_point][i]] = new_ap + temp_f[i] * perp_v

def _get_rot_v(v, n=1, init_pv=None, mrot=210.):
    if init_pv is None:
        temp = abs(v[0])
        pos = 0
        for i in range(1, 3):
            if abs(v[i]) < temp:
                temp = abs(v[i])
                pos = i
        axis = [0., 0., 0.]
        axis[pos] = 1
        pv = rot(v, 90, axis)
        pv[pos] = 0
    else:
        pv = np.array(init_pv)
    pv /= np.linalg.norm(pv)
    rot_vs = np.array(pv)
    if n > 1:
        rot_angles = np.linspace(0, mrot, n)
        for a in rot_angles[1:]:
            rot_vs = np.stack([rot_vs, rot(pv, a, v)])
    return rot_vs

def add_ring_path_3D(path2add, coord_df, anchor_point, temp_v, new_ap, v_fact, **kwargs):
    def _get_pyramid(n):
        if n%2:
            n = n // 2
            a = np.linspace(0,n, n+1)
            a = np.append(a, a[1::-1])
        else:
            n = n // 2
            a = np.linspace(0,n-1, n)
            a = np.append(a, a[::-1])
        return a
    v_fact = kwargs.get('v_fact')
    rnd_fact = v_fact * kwargs.get('rnd_fact', 0.2)
    if len(anchor_point) == 2:
        perp_v = coord_df.loc[list(anchor_point)].diff().loc[anchor_point[1]].values
        perp_v /= np.linalg.norm(perp_v)
    else:
        perp_v = None
    Nat2add = len(path2add[anchor_point])
    rot_vs = _get_rot_v(temp_v, Nat2add, perp_v)
    temp_pyr = _get_pyramid(Nat2add)
    new_positions = np.stack([new_ap]*Nat2add) + temp_pyr * temp_v + rot_vs * v_fact
    new_positions = new_positions + np.random.random(new_positions.shape) * rnd_fact
    for i in range(Nat2add):
        coord_df.loc[path2add[anchor_point][i]] = new_positions[i]

def generate_new_coordinates(path2add, anchor_points, coord_df, v_fact=0.2, **kwargs):
    for anchor_point in anchor_points:
        temp_anch_point = np.average(coord_df.loc[list(anchor_point)], 0)
        temp_v = np.empty((0, coord_df.shape[1]))# make an empty array of an appropriate shape
        for i in range(len(anchor_point)): # anchor_point can be a list of len of 2 (in a ring)
            if anchor_points[anchor_point][i]:
                anchor_p = anchor_point[i]
                temp_N = len(anchor_points[anchor_point][i])
                anchor_p_other = list(anchor_points[anchor_point][i])
                v_diff = coord_df.loc[[anchor_p] * temp_N].values - coord_df.loc[anchor_p_other].values
                temp_v = np.vstack([temp_v, v_diff])
        if not temp_v.any():
            temp_v = [np.ones(coord_df.shape[1])]
        temp_v = np.average(temp_v, 0)
        temp_v = temp_v / np.linalg.norm(temp_v)
        temp_v = temp_v * v_fact
        new_ap = temp_anch_point + temp_v
        if len(path2add[anchor_point]) == 1:
            coord_df.loc[path2add[anchor_point][0]] = new_ap
        else:
            temp_fnc = {1:add_ring_path_simple, 2:add_ring_path_3D}[kwargs.get('fnc_ring_path', 2)]
            temp_fnc(path2add, coord_df, anchor_point, temp_v, new_ap, v_fact, **kwargs)

def get_aligned_coord(v1, v2, v2_align_on=None, weights=None, v1_0=None, cog_1=None):
    """
    :param v1: vector 1 (template for alignment)
    :param v2: vector 2 (vector that is aligned)
    :param v2_align_on: a subsection of v2 used for alignment (same dimension as v1). if None, v2 is used
    :param weights: weights of different atoms for alignment
    :param v1_0: v1 translated to 0,0,0. if None, v1_0 = v1 - cog_1
    :param cog_1: center of geomerty of v1. if None, cog_1 = np.average(v1, axis=0, weights=weights)
    :return:
        v2_aligned
    """
    if cog_1 is None:
        cog_1 = np.average(v1, axis=0, weights=weights)
    if v1_0 is None:
        v1_0 = v1 - cog_1
    if v2_align_on is not None:
        v2_cog = np.average(v2_align_on, axis=0, weights=weights)
        v2_align_on_0 = v2_align_on - v2_cog
        v2_0 = v2 - v2_cog
    else:
        v2_0 = v2 - np.average(v2, axis=0, weights=weights)
        v2_align_on_0 = v2_0
    rot_M, rssd = _align_vectors(v1_0, v2_align_on_0, weights=weights)
    return rot_M.apply(v2_0) + cog_1
