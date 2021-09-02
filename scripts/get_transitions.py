import numpy as np
from SMArt.md.ana.incl import get_lifetime_trans

"""
Author: Drazen Petrov
Notes: script adopted from the original "make_transitions.py" by Jan Perthold
"""

if __name__ == '__main__':
    #------------------------------------------------------
    import argparse

    parser = argparse.ArgumentParser(description="Analyze timeseries of visited states from EDS simulation.")
    parser.add_argument('-n','--numstates', type=int, default="2", required=True, help="Number of EDS states")
    parser.add_argument('-T','--temperature', type=float, default=300, required=False, help="Temperature")
    parser.add_argument('-i','--input', type=str, default="statetser.dat", required=False, help="Filename to print timeseries of visited states.")
    parser.add_argument('-skip', type=int, default=0, required=False, help="skip first n rows from the timeseries file")
    parser.add_argument('-stride', type=int, default=1, required=False, help="skip every Nth row from the timeseries file")
    parser.add_argument('-time', type=float, default=None, nargs=2, required=False, help="t and dt to generate time points")
    parser.add_argument('-include_first_state', default=False, required=False, action='store_true', help="includes first state in the analysis")
    parser.add_argument('-include_last_state', default=False, required=False, action='store_true', help="includes last state in the analysis")
    parser.add_argument('-state_offset', type=int, default=1, required=False, help="adds offset (conuting from 0) to states")
    args=parser.parse_args()
    #------------------------------------------------------
    data_tser = np.loadtxt(args.input, skiprows=args.skip)
    if args.stride!=1:
        data_tser = data_tser[::args.stride, :]
    t, states_tser = data_tser.T
    states_tser = states_tser.astype(int)
    if args.time:
        t = np.arange(args.time[0], args.time[0] + args.time[1] * len(t), args.time[1])

    if args.stride!=1:
        states_tser = states_tser[::args.stride]

    state_life_times, state_trans = get_lifetime_trans(states_tser, args.numstates, time_tser=t, 
        include_first_state=args.include_first_state, include_last_state=args.include_last_state)
    
    ########## output ##########
    print("###")
    print("Average State Lifetimes:")
    for st in range(args.numstates):
        temp_txt = 'State {:}:\t'.format(st + args.state_offset)
        if state_life_times[st]:
            print(temp_txt + '{:.2f}'.format(np.mean(state_life_times[st])))
        else:
            print(temp_txt + 'NaN')
    t_tot = 0
    first_non_zero_state = None
    dGs = []

    print("\n###")
    print("State Sampling Times:")
    for st in range(args.numstates):
        temp_sum = np.sum(state_life_times[st])
        if temp_sum != 0:
            if first_non_zero_state is None:
                first_non_zero_state = st
                first_non_zero_state_t_tot = temp_sum
            temp_dG = -8.314 * args.temperature * np.log(temp_sum / first_non_zero_state_t_tot) / 1000
            dGs.append(temp_dG)
        else:
            dGs.append(None)
        t_tot += temp_sum
        print("State {:}:\t{:.2f}".format(st + args.state_offset, temp_sum))

    print("\n###\nTotal Sampling Time:\t{:.2f}".format(t_tot))

    if first_non_zero_state is not None:
        print("\n###")
        temp_txt = "Free-Energy Differences to State {:} (kJ/mol, at {:} K):"
        print(temp_txt.format(first_non_zero_state + args.state_offset, args.temperature))
        for st in range(args.numstates):
            temp_txt = "State {:}:\t".format(st + args.state_offset)
            if dGs[st] is None:
                print(temp_txt + 'NaN')
            else:
                print(temp_txt + '{:.2f}'.format(dGs[st]))

    print("\n###")
    print("Transitions to:\tFrom:")
    print(state_trans)
