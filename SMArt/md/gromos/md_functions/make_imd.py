_disres_block = """DISTANCERES
#   NTDIR -2..2 controls distance restraining
#         0: no distrance restraining (default)
#         1: instantaneous, using force constant CDIR
#         2: instantaneous, using force constant CDIR x W0
#        -1: time-averaged, using force constant CDIR
#        -2: time-averaged, using force constant CDIR x W0
#  NTDIRA 0,1 controls values for initial distance averages
#         0: generate initial averages
#         1: read from configuration
#    CDIR >= 0.0 force constant for distance restraining
#    DIR0 > 0.0 distance offset in restraining function
#  TAUDIR >= 0.0 coupling time for time averaging
# FORCESCALE 0..2 controls approximation of force scaling
#         0: approximate d<r>/dr = 1
#         1: approximate d<r>/dr = (1.0 - exp(-Dt/tau))
#         2: use d<r>/dr = (1.0 - exp(-Dt/tau))*(<r>/r)^4
#    VDIR 0,1 controls contribution to virial
#         0: no contribution
#         1: distance restraints contribute to virial
#  NTWDIR >= 0 write every NTWDIRth step dist. restr. information to external file
#   NTDIR  NTDIRA    CDIR    DIR0  TAUDIR  FORCESCALE  VDIR  NTWDIR
        2       0     1.0     1.0     1.0           0     0       0
END"""


def _add_imd_blocks(**kwargs):
    s = ''
    if kwargs.get('dsr', False):
        s += _disres_block
    return s + '\n'


def get_emin_missing_at_imd(n_atoms, box_type='0', **kwargs):
    s = """TITLE
Energy minimisation of 1pyo_l1A in vacuum
Created 2017-05-10 13:24:42
END
SYSTEM
#      NPM      NSM
         1      0
END
ENERGYMIN
#     NTEM    NCYC    DELE    DX0     DXM   NMIN   FLIM
         1       1     0.001     0.01    0.05   1000   0.0
END
STEP
#   NSTLIM         T        DT
   5000      0.202      0.002
END
INITIALISE
#    NTIVEL   NTISHK   NTINHT    NTINHB    NTISHI  NTIRTC     NTICOM   NTISTI      IG     TEMPI
      0       0        0         0         1       0          0        0       184930     298.15
END
BOUNDCOND
#      NTB    NDFMIN
       """ + box_type + """      1
END
FORCE
#      NTF array
# bonds    angles    imp.     dihe     charge nonbonded
# H        H         H        H
     1        1         1        1     0  0
# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)
    1 """ + str(n_atoms) + """
END
# use the shake algorithm to constrain the bond lengths.
CONSTRAINT
#      NTC       NTCP   NTCP0(1)     NTCS      NTCS0(1)
         0          1    0.00010        1      0.00010
END
PAIRLIST
#    algorithm: standard (0) (gromos96 like pairlist)
#                     grid (1) (XX grid pairlist)
#    SIZE:       grid cell size (or auto = 0.5 * RCUTP)
#    TYPE:       chargegoup (0) (chargegroup based cutoff)
#                     atomic (1) (atom based cutoff)
#
#    algorithm    NSNB    RCUTP    RCUTL    SIZE    TYPE
            0      5      0.8      1.4     0.4     0
END
NONBONDED
# NLRELE    APPAK      RCRF     EPSRF  NSLFEXCL
       1      0.0       1.4          1     1
# NSHAPE   ASHAPE    NA2CLC   TOLA2   EPSLS
       -1       1.4        2   1e-10       0
# NKX    NKY   NKZ    KCUT
   10     10    10     100
# NGX   NGY   NGZ  NASORD  NFDORD   NALIAS  NSPORD
   32    32    32       3       2        3       4
# NQEVAL   FACCUR   NRDGRD   NWRGRD   NLRLJ    SLVDNS
  100000      1.6        0        0       0      33.3
END
COMTRANSROT
#   NSCM
    0
END
PRINTOUT
#NTPR: print out energies, etc. every NTPR steps
#NTPP: =1 perform dihedral angle transition monitoring
#     NTPR          NTPP
        10          0
END
POSITIONRES
# values for NTPOR:
#  0: no position re(con)straining
#  1: use CPOR
#  2: use SPOR/ ATOMIC B-FACTORS
#  3: position constraining
# NTPOR       NTPORB     NTPORS      CPOR
      3            1          0     2.5E4
END
COVALENTFORM
#    NTBBH    NTBAH     NTBDN
         0         0         0
END
WRITETRAJ
#    NTWX     NTWSE      NTWV      NTWF      NTWE      NTWG      NTWB
      0         0         0         0         0         0         0
END
INNERLOOP
        0  0  0
END
ROTTRANS
#      RTC   RTCLAST
       0      """ + str(n_atoms) + """
END
"""
    s += _add_imd_blocks(**kwargs)
    return s


def get_emin_vac_imd(n_atoms, **kwargs):
    s = """TITLE
Energy minimisation of 1pyo_l1A in vacuum
Created 2017-05-10 13:24:42
END
SYSTEM
#      NPM      NSM
         1      0
END
ENERGYMIN
#     NTEM    NCYC    DELE    DX0     DXM   NMIN   FLIM
         1       1     0.001     0.01    0.05   1000   0.0
END
STEP
#   NSTLIM         T        DT
   5000         0.0      0.002
END
INITIALISE
#    NTIVEL   NTISHK   NTINHT    NTINHB    NTISHI  NTIRTC     NTICOM   NTISTI      IG     TEMPI
      0       0        0         0         1       0          0        0       184930     298.15
END
BOUNDCOND
#      NTB    NDFMIN
       0      1
END
FORCE
#      NTF array
# bonds    angles    imp.     dihe     charge nonbonded
# H        H         H        H
     0        1         1        1     1  1
# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)
    1 """ + str(n_atoms) + """
END
# use the shake algorithm to constrain the bond lengths.
CONSTRAINT
#      NTC       NTCP   NTCP0(1)     NTCS      NTCS0(1)
         3          1    0.00010        1      0.00010
END
PAIRLIST
#    algorithm: standard (0) (gromos96 like pairlist)
#                     grid (1) (XX grid pairlist)
#    SIZE:       grid cell size (or auto = 0.5 * RCUTP)
#    TYPE:       chargegoup (0) (chargegroup based cutoff)
#                     atomic (1) (atom based cutoff)
#
#    algorithm    NSNB    RCUTP    RCUTL    SIZE    TYPE
            0      5      0.8      1.4     0.4     0
END
NONBONDED
# NLRELE    APPAK      RCRF     EPSRF  NSLFEXCL
       1      0.0       1.4          1     1
# NSHAPE   ASHAPE    NA2CLC   TOLA2   EPSLS
       -1       1.4        2   1e-10       0
# NKX    NKY   NKZ    KCUT
   10     10    10     100
# NGX   NGY   NGZ  NASORD  NFDORD   NALIAS  NSPORD
   32    32    32       3       2        3       4
# NQEVAL   FACCUR   NRDGRD   NWRGRD   NLRLJ    SLVDNS
  100000      1.6        0        0       0      33.3
END
COMTRANSROT
#   NSCM
    0
END
PRINTOUT
#NTPR: print out energies, etc. every NTPR steps
#NTPP: =1 perform dihedral angle transition monitoring
#     NTPR          NTPP
        10          0
END
POSITIONRES
# values for NTPOR:
#  0: no position re(con)straining
#  1: use CPOR
#  2: use SPOR/ ATOMIC B-FACTORS
#  3: position constraining
# NTPOR       NTPORB     NTPORS      CPOR
      0            1          0     2.5E4
END
COVALENTFORM
#    NTBBH    NTBAH     NTBDN
         0         0         0
END
WRITETRAJ
#    NTWX     NTWSE      NTWV      NTWF      NTWE      NTWG      NTWB
      0         0         0         0         0         0         0
END
INNERLOOP
        0  0  0
END
ROTTRANS
#      RTC   RTCLAST
       0      """ + str(n_atoms) + """
END
"""
    s += _add_imd_blocks(**kwargs)
    return s


def get_emin_water_imd(n_atoms, nW, **kwargs):
    s = """TITLE
steepest descent energy minimization of the peptide in water
END
ENERGYMIN
#     NTEM    NCYC    DELE    DX0     DXM   NMIN   FLIM
         1       1      0.1   0.01    0.05   1000    0.0
END
SYSTEM
#      NPM      NSM
         1      """ + str(nW) + """
END
INITIALISE
#    NTIVEL   NTISHK   NTINHT    NTINHB    NTISHI  NTIRTC     NTICOM   NTISTI      IG     TEMPI
         0         0        0         0         1       0          0        0  210185     300.0
END
# do 2000 steps
STEP
#   NSTLIM         T        DT
      1000       0.0     0.002
END
# do it with rectangular periodic boundary conditions
BOUNDCOND
#      NTB    NDFMIN
         1         1
END
# every 10 steps print the energy in the output file.
PRINTOUT
#NTPR: print out energies, etc. every NTPR steps
#NTPP: =1 perform dihedral angle transition monitoring
#     NTPR          NTPP
        10             0
END
# use the shake algorithm to constrain the bond lengths.
CONSTRAINT
#      NTC       NTCP   NTCP0(1)     NTCS      NTCS0(1)
         3          1    0.00010        1      0.00010
END
FORCE
#      NTF array
# bonds    angles    imp.     dihe     charge nonbonded
# H        H         H        H
     0        1         1        1     1  1
# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)
     2      """ + str(n_atoms) + """      """ + str(nW * 3 + n_atoms) + """
END
# with rectangular periodic boundary conditions we may use
# the grid based pairlist generation
PAIRLIST
#	algorithm: standard (0) (gromos96 like pairlist)
#		             grid (1) (XX grid pairlist)
#	SIZE:	   grid cell size (or auto = 0.5 * RCUTP)
#	TYPE:	   chargegoup (0) (chargegroup based cutoff)
#			         atomic (1) (atom based cutoff)
#
#	algorithm	NSNB	RCUTP	RCUTL	SIZE	TYPE
	0               5	  0.8	  1.4	0.4     0
END
NONBONDED
# NLRELE    APPAK      RCRF     EPSRF  NSLFEXCL
       1      0.0       1.4        61     1
# NSHAPE   ASHAPE    NA2CLC   TOLA2   EPSLS
      -1       1.4        2   1e-10       0
# NKX    NKY   NKZ    KCUT
   10     10    10     100
# NGX   NGY   NGZ  NASORD  NFDORD   NALIAS  NSPORD
   32    32    32       3       2        3       4
# NQEVAL   FACCUR   NRDGRD   NWRGRD   NLRLJ    SLVDNS
  100000      1.6        0        0       0      33.3
END
POSITIONRES
# values for NTPOR:
#  0: no position re(con)straining
#  1: use CPOR
#  2: use SPOR/ ATOMIC B-FACTORS
#  3: position constraining
# NTPOR       NTPORB     NTPORS      CPOR
      0            1          0     2.5E4
END
"""
    s += _add_imd_blocks(**kwargs)
    return s


def get_eq_imd(n_atoms, nW, **kwargs):
    gpu_flag = kwargs.get('gpu_flag', True)
    if gpu_flag:
        inner_loop = '         4         0         1'
    else:
        inner_loop = '         2         1         0'
    NSCM = kwargs.get('NSCM', 1000)
    s = """TITLE
equilibration of the peptide in water
END
# we have 1 solute and 908 solvent molecules
SYSTEM
#      NPM      NSM
         1      """ + str(nW) + """
END
# most of this block is overwritten by mkscript.
INITIALISE
#    NTIVEL   NTISHK   NTINHT    NTINHB    NTISHI  NTIRTC     NTICOM   NTISTI      IG     TEMPI
         0         0        0         0         1       0          0        0  210185       0.0
END
# do 50000 steps
STEP
#   NSTLIM         T        DT
     50000       0.0     0.002
END
# do it with rectangular periodic boundary conditions
BOUNDCOND
#      NTB     NDFMIN
         1         3
END
# couple the temperature, the temperatures are overwritten by mkscript.
MULTIBATH
# ALGORITHM:
#      weak-coupling(0):      use weak-coupling scheme
#      nose-hoover(1):        use Nose Hoover scheme
#      nose-hoover-chains(2): use Nose Hoover chains scheme
# NUM: number of chains in Nose Hoover chains scheme
#      !! only specify NUM when needed !!
# NBATHS: number of temperature baths to couple to
#          ALGORITHM
                   0
#  NBATHS
         2
# TEMP0(1 ... NBATHS)  TAU(1 ... NBATHS)
       300     0.1     300     0.1
#   DOFSET: number of distiguishable sets of d.o.f.
         2
# LAST(1 ... DOFSET)  COMBATH(1 ... DOFSET)  IRBATH(1 ... DOFSET)
       """ + str(n_atoms) + """         1         1       """ + str(nW * 3 + n_atoms) + """         2         2
END
PRESSURESCALE
# COUPLE   SCALE    COMP    TAUP  VIRIAL
       2       1 0.0004575      0.5       2
# SEMI
       1       1       1
# PRES0(1...3,1...3)
 0.06102       0       0
       0 0.06102       0
       0       0 0.06102
END
# every 1000 step we remove only the translational com motion
COMTRANSROT
#   NSCM
    """ + str(NSCM) + """
END
ROTTRANS
#      RTC   RTCLAST
       0      """ + str(n_atoms) + """
END
INNERLOOP
""" + inner_loop + """
END
COVALENTFORM
# NTBBH: 0,1 controls bond-stretching potential
#        0: quartic form (default)
#        1: harmonic form
# NTBAH: 0,1 controls bond-angle bending potential
#        0: cosine-harmonic (default)
#        1: harmonic
# NTBDN: 0,1 controls torsional dihedral potential
#        0: arbitrary phase shifts (default)
#        1: phase shifts limited to 0 and 180 degrees.
#   NTBBH    NTBAH    NTBDN
        0        0        0
END
# every 250 steps write the energy and coordinates to the
# trajectory
WRITETRAJ
# NTWSE = configuration selection parameter
# =0: write normal trajectory
# >0: chose min energy for writing configurations
#     NTWX     NTWSE      NTWV      NTWF    NTWE      NTWG      NTWB
      5000         0         0         0    5000         0         0
END
# every 250 steps print the energy in the output file.
PRINTOUT
#NTPR: print out energies, etc. every NTPR steps
#NTPP: =1 perform dihedral angle transition monitoring
#     NTPR     NTPP
      5000        0
END
# calculate the energies between the peptide, the ions and the solvent.
FORCE
 # NTF(1..6): 0,1 determines terms used in force calculation
#             0: do not include terms
#             1: include terms
# NEGR: ABS(NEGR): number of energy groups
#             > 0: use energy groups
#             < 0: use energy and force groups
# NRE(1..NEGR): >= 1.0 last atom in each energy group
# NTF(1) NTF(2) NTF(3) NTF(4) NTF(5)        NTF(6)
# bonds     angles    improper  dihedral  electrostatic vdW
  0         1         1         1         1             1
# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)
     2      """ + str(n_atoms) + """    """ + str(nW * 3 + n_atoms) + """
END
# use the shake algorithm to constrain the bond lengths.
CONSTRAINT
#      NTC       NTCP   NTCP0(1)     NTCS      NTCS0(1)
         3          1    0.00010       4
END
# use grid based pairlist generation to speed up
PAIRLIST
#	algorithm: standard(0) (gromos96 like pairlist)
#		             grid(1) (XX grid pairlist)
#	SIZE:	   grid cell size (or auto = 0.5 * RCUTP)
#	TYPE:	   chargegoup(0) (chargegroup based cutoff)
#			         atomic(1) (atom based cutoff)
#
#	algorithm	  NSNB	RCUTP	RCUTL	  SIZE	TYPE
	        1	   5	  0.8	  1.4	   0.4	   0
END
# Longrange reaction field correction
NONBONDED
# NLRELE    APPAK      RCRF     EPSRF   NSLFEXCL
       1      0.0       1.4        61    1
# NSHAPE   ASHAPE    NA2CLC   TOLA2   EPSLS
      -1       1.4        2   1e-10       0
# NKX    NKY   NKZ    KCUT
   10     10    10     100
# NGX   NGY   NGZ  NASORD  NFDORD   NALIAS  NSPORD
   32    32    32       3       2        3       4
# NQEVAL   FACCUR   NRDGRD   NWRGRD   NLRLJ    SLVDNS
  100000      1.6        0        0       0      33.3
END
POSITIONRES
# values for NTPOR:
#  0: no position re(con)straining
#  1: use CPOR
#  2: use SPOR/ ATOMIC B-FACTORS
#  3: position constraining
# NTPOR       NTPORB     NTPORS      CPOR
      1            1          0     2.5E4
END
"""
    s += _add_imd_blocks(**kwargs)
    return s


def get_ext_TI_imd(n_atoms, nW, bond_const=False, **kwargs):
    if bond_const:
        const_txt = """#      NTC       NTCP   NTCP0(1)     NTCS      NTCS0(1)
         3          1    0.00010       4"""
        force_txt = '0'
    else:
        const_txt = """#      NTC       NTCP   NTCP0(1)     NTCS      NTCS0(1)
         4          1    0.00010       4"""
        force_txt = '1'
    s = """TITLE
extended TI
END
# we have 1 solute and 908 solvent molecules
SYSTEM
#      NPM      NSM
         1      """ + str(nW) + """
END
# most of this block is overwritten by mkscript.
INITIALISE
#    NTIVEL   NTISHK   NTINHT    NTINHB    NTISHI  NTIRTC     NTICOM   NTISTI      IG     TEMPI
         0         0        0         0         1       0          0        0  210185       0.0
END
# do 50000 steps
STEP
#   NSTLIM         T        DT
     50000       0.0     0.002
END
# do it with rectangular periodic boundary conditions
BOUNDCOND
#      NTB     NDFMIN
         1         3
END
# couple the temperature, the temperatures are overwritten by mkscript.
MULTIBATH
# ALGORITHM:
#      weak-coupling(0):      use weak-coupling scheme
#      nose-hoover(1):        use Nose Hoover scheme
#      nose-hoover-chains(2): use Nose Hoover chains scheme
# NUM: number of chains in Nose Hoover chains scheme
#      !! only specify NUM when needed !!
# NBATHS: number of temperature baths to couple to
#          ALGORITHM
                   0
#  NBATHS
         2
# TEMP0(1 ... NBATHS)  TAU(1 ... NBATHS)
       300     0.1     300     0.1
#   DOFSET: number of distiguishable sets of d.o.f.
         2
# LAST(1 ... DOFSET)  COMBATH(1 ... DOFSET)  IRBATH(1 ... DOFSET)
       """ + str(n_atoms) + """         1         1       """ + str(nW * 3 + n_atoms) + """         2         2
END
PRESSURESCALE
# COUPLE   SCALE    COMP    TAUP  VIRIAL
       2       1 0.0004575      0.5       2
# SEMI
       1       1       1
# PRES0(1...3,1...3)
 0.06102       0       0
       0 0.06102       0
       0       0 0.06102
END
# every 1000 step we remove only the translational com motion
COMTRANSROT
#   NSCM
    1000
END
COVALENTFORM
# NTBBH: 0,1 controls bond-stretching potential
#        0: quartic form (default)
#        1: harmonic form
# NTBAH: 0,1 controls bond-angle bending potential
#        0: cosine-harmonic (default)
#        1: harmonic
# NTBDN: 0,1 controls torsional dihedral potential
#        0: arbitrary phase shifts (default)
#        1: phase shifts limited to 0 and 180 degrees.
#   NTBBH    NTBAH    NTBDN
        0        0        0
END
FORCE
 # NTF(1..6): 0,1 determines terms used in force calculation
#             0: do not include terms
#             1: include terms
# NEGR: ABS(NEGR): number of energy groups
#             > 0: use energy groups
#             < 0: use energy and force groups
# NRE(1..NEGR): >= 1.0 last atom in each energy group
# NTF(1) NTF(2) NTF(3) NTF(4) NTF(5)        NTF(6)
# bonds     angles    improper  dihedral  electrostatic vdW
  """ + force_txt + """         1         1         1         1             1
# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)
     2      """ + str(n_atoms) + """    """ + str(nW * 3 + n_atoms) + """
END
# use the shake algorithm to constrain the bond lengths.
CONSTRAINT
""" + const_txt + """
END
INNERLOOP
         2         1         0
END
# use grid based pairlist generation to speed up
PAIRLIST
#	algorithm: standard(0) (gromos96 like pairlist)
#		             grid(1) (XX grid pairlist)
#	SIZE:	   grid cell size (or auto = 0.5 * RCUTP)
#	TYPE:	   chargegoup(0) (chargegroup based cutoff)
#			         atomic(1) (atom based cutoff)
#
#	algorithm	  NSNB	RCUTP	RCUTL	  SIZE	TYPE
	        1	   5	  0.8	  1.4	   0.4	   0
END
# Longrange reaction field correction
NONBONDED
# NLRELE    APPAK      RCRF     EPSRF   NSLFEXCL
       1      0.0       1.4        61    1
# NSHAPE   ASHAPE    NA2CLC   TOLA2   EPSLS
       3       1.4        2   1e-10       0
# NKX    NKY   NKZ    KCUT
   10     10    10     100
# NGX   NGY   NGZ  NASORD  NFDORD   NALIAS  NSPORD
   32    32    32       3       2        3       4
# NQEVAL   FACCUR   NRDGRD   NWRGRD   NLRLJ    SLVDNS
  100000      1.6        0        0       0      33.3
END
PRINTOUT
#NTPR: print out energies, etc. every NTPR steps
#NTPP: =1 perform dihedral angle transition monitoring
#     NTPR      NTPP
      1000         0
END
WRITETRAJ
#    NTWX     NTWSE      NTWV      NTWF      NTWE      NTWG      NTWB
     1000         0         0         0        20        20         0
END
PERTURBATION
#    NTG: 0..1 controls use of free-energy calculation.
#         0: no free-energy calculation (default)
#         1: calculate dH/dRLAM
#  NRDGL: 0,1 controls reading of initial value for RLAM.
#         0: use initial RLAM parameter from PERTURBATION block
#         1: read from configuration
#   RLAM: 0.0..1.0 initial value for lambda
#  DLAMT: >= 0.0 rate of lambda increase in time.
# ALPHLJ: >= 0.0 Lennard-Jones soft-core parameter
#  ALPHC: >= 0.0 Coulomb-RF soft-core parameter
#   NLAM: > 0 power dependence of lambda coupling
# NSCALE: 0..2 controls use of interaction scaling
#         0: no interaction scaling
#         1: interaction scaling
#         2: perturbation for all atom pairs with scaled
#            interactions. No perturbation for others.
#
#     NTG   NRDGL    RLAM   DLAMT
        1       0     0.0     0.0
#  ALPHLJ   ALPHC    NLAM  NSCALE
      0.5     0.5       1       0
END
PRECALCLAM
# NRLAM   0  : off
#         >1 : precalculating energies for NRLAM extra lambda values
# MINLAM  between 0 and 1: minimum lambda value to precalculate energies
# MAXLAM  between MINLAM and 1: maximum lambda value to precalculate energies
# NRLAM   MINLAM   MAXLAM
     81      0.0      1.0
END
"""
    s += _add_imd_blocks(**kwargs)
    return s


def get_short_TI_imd(n_atoms, nW, **kwargs):
    s = """TITLE
equilibration of the peptide in water
END
# we have 1 solute and 908 solvent molecules
SYSTEM
#      NPM      NSM
         1      """ + str(nW) + """
END
# most of this block is overwritten by mkscript.
INITIALISE
#    NTIVEL   NTISHK   NTINHT    NTINHB    NTISHI  NTIRTC     NTICOM   NTISTI      IG     TEMPI
         0         0        0         0         1       0          0        0  210185       0.0
END
# do 50000 steps
STEP
#   NSTLIM         T        DT
      5000       0.0     0.002
END
# do it with rectangular periodic boundary conditions
BOUNDCOND
#      NTB     NDFMIN
         1         3
END
# couple the temperature, the temperatures are overwritten by mkscript.
MULTIBATH
# ALGORITHM:
#      weak-coupling(0):      use weak-coupling scheme
#      nose-hoover(1):        use Nose Hoover scheme
#      nose-hoover-chains(2): use Nose Hoover chains scheme
# NUM: number of chains in Nose Hoover chains scheme
#      !! only specify NUM when needed !!
# NBATHS: number of temperature baths to couple to
#          ALGORITHM
                   0
#  NBATHS
         2
# TEMP0(1 ... NBATHS)  TAU(1 ... NBATHS)
       300     0.1     300     0.1
#   DOFSET: number of distiguishable sets of d.o.f.
         2
# LAST(1 ... DOFSET)  COMBATH(1 ... DOFSET)  IRBATH(1 ... DOFSET)
       """ + str(n_atoms) + """         1         1       """ + str(nW * 3 + n_atoms) + """         2         2
END
PRESSURESCALE
# COUPLE   SCALE    COMP    TAUP  VIRIAL
       2       1 0.0004575      0.5       2
# SEMI
       1       1       1
# PRES0(1...3,1...3)
 0.06102       0       0
       0 0.06102       0
       0       0 0.06102
END
# every 1000 step we remove only the translational com motion
COMTRANSROT
#   NSCM
    1000
END
COVALENTFORM
# NTBBH: 0,1 controls bond-stretching potential
#        0: quartic form (default)
#        1: harmonic form
# NTBAH: 0,1 controls bond-angle bending potential
#        0: cosine-harmonic (default)
#        1: harmonic
# NTBDN: 0,1 controls torsional dihedral potential
#        0: arbitrary phase shifts (default)
#        1: phase shifts limited to 0 and 180 degrees.
#   NTBBH    NTBAH    NTBDN
        0        0        0
END
# every 250 steps write the energy and coordinates to the
# trajectory
WRITETRAJ
# NTWSE = configuration selection parameter
# =0: write normal trajectory
# >0: chose min energy for writing configurations
#     NTWX     NTWSE      NTWV      NTWF    NTWE      NTWG      NTWB
      5000         0         0         0    5000         0         0
END
# every 250 steps print the energy in the output file.
PRINTOUT
#NTPR: print out energies, etc. every NTPR steps
#NTPP: =1 perform dihedral angle transition monitoring
#     NTPR     NTPP
      5000        0
END
# calculate the energies between the peptide, the ions and the solvent.
FORCE
 # NTF(1..6): 0,1 determines terms used in force calculation
#             0: do not include terms
#             1: include terms
# NEGR: ABS(NEGR): number of energy groups
#             > 0: use energy groups
#             < 0: use energy and force groups
# NRE(1..NEGR): >= 1.0 last atom in each energy group
# NTF(1) NTF(2) NTF(3) NTF(4) NTF(5)        NTF(6)
# bonds     angles    improper  dihedral  electrostatic vdW
  0         1         1         1         1             1
# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)
     2      """ + str(n_atoms) + """    """ + str(nW * 3 + n_atoms) + """
END
# use the shake algorithm to constrain the bond lengths.
CONSTRAINT
#      NTC       NTCP   NTCP0(1)     NTCS      NTCS0(1)
         3          1    0.00010       4
END
# use grid based pairlist generation to speed up
PAIRLIST
#	algorithm: standard(0) (gromos96 like pairlist)
#		             grid(1) (XX grid pairlist)
#	SIZE:	   grid cell size (or auto = 0.5 * RCUTP)
#	TYPE:	   chargegoup(0) (chargegroup based cutoff)
#			         atomic(1) (atom based cutoff)
#
#	algorithm	  NSNB	RCUTP	RCUTL	  SIZE	TYPE
	        1	   5	  0.8	  1.4	   0.4	   0
END
# Longrange reaction field correction
NONBONDED
# NLRELE    APPAK      RCRF     EPSRF   NSLFEXCL
       1      0.0       1.4        61    1
# NSHAPE   ASHAPE    NA2CLC   TOLA2   EPSLS
      -1       1.4        2   1e-10       0
# NKX    NKY   NKZ    KCUT
   10     10    10     100
# NGX   NGY   NGZ  NASORD  NFDORD   NALIAS  NSPORD
   32    32    32       3       2        3       4
# NQEVAL   FACCUR   NRDGRD   NWRGRD   NLRLJ    SLVDNS
  100000      1.6        0        0       0      33.3
END
INNERLOOP
         2         1         0
END
POSITIONRES
# values for NTPOR:
#  0: no position re(con)straining
#  1: use CPOR
#  2: use SPOR/ ATOMIC B-FACTORS
#  3: position constraining
# NTPOR       NTPORB     NTPORS      CPOR
      1            1          0        50
END
PERTURBATION
#    NTG: 0..1 controls use of free-energy calculation.
#         0: no free-energy calculation (default)
#         1: calculate dH/dRLAM
#  NRDGL: 0,1 controls reading of initial value for RLAM.
#         0: use initial RLAM parameter from PERTURBATION block
#         1: read from configuration
#   RLAM: 0.0..1.0 initial value for lambda
#  DLAMT: >= 0.0 rate of lambda increase in time.
# ALPHLJ: >= 0.0 Lennard-Jones soft-core parameter
#  ALPHC: >= 0.0 Coulomb-RF soft-core parameter
#   NLAM: > 0 power dependence of lambda coupling
# NSCALE: 0..2 controls use of interaction scaling
#         0: no interaction scaling
#         1: interaction scaling
#         2: perturbation for all atom pairs with scaled
#            interactions. No perturbation for others.
#
#     NTG   NRDGL    RLAM   DLAMT
        1       0     0.0     0.1
#  ALPHLJ   ALPHC    NLAM  NSCALE
      0.5     0.5       1       0
END
"""
    s += _add_imd_blocks(**kwargs)
    return s
