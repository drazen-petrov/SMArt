from SMArt.md.incl import _gr_title_gm_system
_maketopv = '_gr_maketopv'
_topv = '_gr_topv'
_physconst = '_gr_physconst'
nexcl = 'nexcl'

_gromos_block_names = dict(FORCEFIELD='ff_name', MAKETOPVERSION=_maketopv, TOPVERSION=_topv,
                           PHYSICALCONSTANTS=_physconst, LINKEXCLUSIONS=nexcl)

_physconst_txt = """# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)
138.9354
# HBAR: Planck's constant HBAR = H/(2* PI)
0.0635078
# SPDL: Speed of light (nm/ps)
299792.458
# BOLTZ: Boltzmann's constant kB
0.00831441"""

_gromos_dvalues = {'ff_name': 'GROMOS_ff', _maketopv: '1.0', _topv: '2.0', nexcl: '2', _physconst: _physconst_txt}
