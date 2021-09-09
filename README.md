Simulation and Modelling Art (SMArt)
=====
[![Documentation](https://img.shields.io/badge/Documentation-here-white.svg)](https://drazen-petrov.github.io/SMArt/)


This toolkit is aimed at facilitating setup and analysis of MD simulations, e.g. perturbation free energy calculations.

## Installation
Just clone this repository and you are ready to go

## Usage
At the moment, the package supports GROMOS and GROMACS file formats

### I/O functionality
parsing

```python
from SMArt.md import parse_top, parse_ff
ff_gr = parse_ff(gromos_ifp_file) # by default, format_type = 'gr'
ff_gm = parse_ff (gromacs_ff_itp_file, format_type = 'gm')

top_gr = parse_top(gromos_top_file) # by default, format_type = 'gr'
top_gm = parse_top(gromacs_top_file, format_type = 'gm')
```
writing
```python
ff_gr.write_ff(output_file)
ff_gm.write_ff(output_file, format_type = 'gm')

top_gr.write_top(output_file)
top_gm.write_top(output_file, format_type = 'gm') # this will produce a single topology file containing FF information
gm_top_kwargs = dict(format_type='gm', flag_if=True, flag_include=True, sep_ff2itp='top_ff.itp', sep_mol2itp=True)
top_gm.write_top(output_file, flag_segment_order=1, flag_defines_first_include=1, **gm_top_kwargs)
# this separates the FF and molecule type information into separate *itp files that are included in the topology
```

### GROMOS <-> GROMACS conversion
use `gr2gm()` function to converte GROMOS to GROMACS format and `gm2gr()` other way around
```python
ff_gr.gr2gm()
top_gr.gr2gm()
ff_gr.write_ff(converted_output_file, format_type='gm')
top_gr.write_top(converted_output_file, format_type='gm')

ff_gm.gm2gr()
top_gm.gm2gr()
ff_gm.write_ff(converted_output_file, format_type='gr')
top_gm.write_top(converted_output_file, format_type='gr')
```

### Perturbation topology preparation
Perturbation topology is based on the maximum common substructure search. This example applies for both point mutation, as well as topologies of 2 different small molecules (e.g. ligands)
```python
from SMArt.md import parse_top
from SMArt import alchemy

format_type = 'gr' # or 'gm'
top_wt = parse_top(wt_top_file, format_type=format_type)
top_mut = parse_top(mut_top_file, format_type=format_type)
# GROMACS
mt_wt, mt_mut = top_gm.molecule_types[mol_id_wt], top_gm_mut.molecule_types[mol_id_mut]
sol = alchemy.point_mutation(mt_wt, mt_mut, ff_dumm = top_gm.get_DUM_type)
sol.toptp.write_itp('toptp.itp', flag_generate_excl_pairs = True)
# GROMOS
sol = alchemy.point_mutation(top_wt, top_mut)
sol.toptp.write_top('toptp.top')
sol.toptp.write_ptp('toptp.ptp')
```

### Simulation update (perturbation free energy calculations)
Automated simulation update scheme - choice of lambda points and simulation time
```python
from SMArt.md.ana import pert_FE
from SMArt.md.gromos.io.ana import read_bar_dhdl, read_exTI
from SMArt.md.gromacs.io.ana import read_bar_data, read_xvg_data

bar_data = {}
for f_path in FE_data_files:
    data, sim_l = read_bar_data(f_path) # gromacs
    data, sim_l = read_bar_dhdl(f_path) # gromos
    pert_FE.combine_bar_dhdl(bar_data, data, sim_l)

initial_t_sim = 1. # initial simulation time per lambda point (1 ns)
LPs_times = dict((temp_lp, initial_t_sim) for temp_lp in bar_data)

seg_score_flag, converged_segments, new_LPs_weights, seg_data_dG_err = pert_FE.update_LPs_times(bar_data)
iter_sim_time = 11. # total simulation time per iteration (if scheme II is used - 11 ns in this case)
new_LPs_times = pert_FE.get_LPs_times(new_LPs_weights, LPs_times, iter_sim_time)
```

More information is provided in `doc` directory of the package<br>
Additional useful scripts are can be found in `scripts` directory of the package

## License
`SMArt` is published under the GNU General Public License [(GPLv3)](https://www.gnu.org/licenses/gpl-3.0.html)

## How to cite?
Please use this reference for citation:<br>
D Petrov. Perturbation Free-Energy Toolkit: An Automated Alchemical Topology Builder. J. Chem. Inf. Model. 2021. ([DOI)](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00428)

