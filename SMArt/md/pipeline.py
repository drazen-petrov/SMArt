from SMArt.incl import os, np, Defaults
from SMArt.md.ana.incl import get_lset, fix_float, Real

def get_slurm_job_settings(**slurm_kw):
    """
    :param slurm_kw: key value dictionary that will add the following lines: '#SBATCH -{:} {:}\n'
    :return: header string
    """
    s = '#!/bin/bash\n'
    for item in slurm_kw.items():
        s+= '#SBATCH -{:} {:}\n'.format(*item)
    return s + '\n'

def write_comm_in_f(f, comm):
    if not comm.endswith('\n'):
        comm+='\n'
    f.write(comm)


def gen_slurm_job(f_path, commands, slurm_kw=None, **kwargs):
    """
    :param f_path: run file
    :param commands: commands to be executed in the run file
    :param slurm_kw: additional slurm settings (dict) - see  get_slurm_job_settings function
    :param kwargs:
        temp_fd: use temp folder to run the job - e.g. /scratch/${SLURM_JOBID}/
        job_fd: job/simulation folder
        flag_cp_all: if temp_job_fd given, this flag ensures that all the data is c/p to/from temp folder - it can also be (False, True) - copies only back to JOBDIR
        pre_commands: set of commands to run before the commands
        post_commands: set of commands to run after the commands
    :return: None
    """
    if slurm_kw is None:
        slurm_kw={}
    f = open(f_path, 'w')
    f.write(get_slurm_job_settings(**slurm_kw))
    temp_fd = kwargs.get('temp_fd')
    job_fd = kwargs.get('job_fd')
    if temp_fd:
        assert job_fd
        f.write('# create temp workdir\n')
        f.write('JOBDIR=' + job_fd + '\n')
        f.write('WORKDIR=' + temp_fd + '\n')
        f.write('mkdir -p ${WORKDIR}\n')
        f.write('cd       ${WORKDIR}\n')
        flag_cp_all = kwargs.get('flag_cp_all')
        if flag_cp_all == True:
            flag_cp_all = (True, True)
        if flag_cp_all == False:
            flag_cp_all = (False, False)
        if flag_cp_all[0]:
            f.write('cp -p ${JOBDIR}/* .\n')
        f.write('\n')
    pre_commands = kwargs.get('pre_commands')
    if pre_commands:
        f.write('# pre-commands\n')
        for comm in pre_commands:
            write_comm_in_f(f, comm)
        f.write('\n')
    f.write('# job commands\n')
    for comm in commands:
        write_comm_in_f(f, comm)
    f.write('\n')
    if temp_fd and flag_cp_all[1]:
        for comm in get_cp_back_commands():
            write_comm_in_f(f, comm)
    post_commands = kwargs.get('post_commands')
    if post_commands:
        f.write('# post-commands\n')
        for comm in post_commands:
            write_comm_in_f(f, comm)
    f.close()
    os.system('chmod a+x '+ f_path)

def get_cp_back_commands():
    commands = []
    commands.append('# copy back from temp workdir\n')
    commands.append('cd ${JOBDIR}\n')
    commands.append('OK=1\n')
    commands.append('cp -r ${WORKDIR}/* ${JOBDIR}/. || OK=0\n')
    commands.append('if `test ${OK} -eq 0`; then\n')
    commands.append('    rm -r ${WORKDIR}\nfi\n\n')
    return commands

def prep_pdb(pdb):
    pass

def eq(top, conf, md_in, eq_jobs):
    pass


class SimulationSet(Defaults):
    @staticmethod
    def split_single_sim(t, t_max):
        assert t == round(t,1)
        if t_max < t:
            N_sims = int(t // t_max)
            if t % t_max > 0.1:N_sims+=1
            t_sim = round(t / N_sims, 1)
            sims = [t_sim] * (N_sims - 1)
            sims.append(round(t - (N_sims - 1) * t_sim, 1))
        else:
            sims = [t]
        return sims

    @staticmethod
    def get_sim_count(sims_count, sim_key, **kwargs):
        """
        based on sims_count (dict), gets the counter for a given simulation type (sim_key)
        also, increases the counter by 1
        :param sims_count: dict (stores the counter)
        :param sim_key: simulation key (associated with the counter)
        :param kwargs:
            flag_increase_count: to increase count by 1 (not done if the commands are just appended in the same job file)
        """
        if sim_key not in sims_count:
            sims_count[sim_key] = 0
        if kwargs.get('flag_increase_count', True):
            sims_count[sim_key] += 1
        return sims_count[sim_key]

    def __fnc2be_implemented(self):
        self.get_sim_file_fd = None # should based on the input simulation, generates a folder; job file name; and simulation name

    def process_sim(self, commands, sim, sim_fd_root, fnc2process, **kwargs):
        #job_info = fnc2process
        job_f, sim_fd, name = self.get_sim_file_fd(sim, sim_fd_root, **kwargs)
        commands.extend(fnc2process(sim, self.sim_set, job_f, sim_fd, name, **kwargs))
        commands.append('\n')
        if sim.get('flag_cp_before_new_sub'):
            commands.extend(get_cp_back_commands())
        # process all jobs to sub
        for new_sim_i in sim.get('sub'):
            new_sim = self.sim_set[new_sim_i]
            new_commands = []
            new_sim_info = self.process_sim(new_commands, new_sim, sim_fd_root, fnc2process, **kwargs)
            slurm_job_kwargs = new_sim.get('slurm_job_kwargs', kwargs)
            gen_slurm_job(new_sim_info[0], new_commands, job_fd=new_sim_info[1], **slurm_job_kwargs) # job_f, commands (kwargs, e.g. slurm_kw)
            temp_job_f_name = os.path.split(new_sim_info[0])[1]
            #if new_sim.get('flag_cp_before_new_sub'):
            #    commands.extend(get_cp_back_commands())
            comm = 'cd ' + new_sim_info[1] + '\n{:} '.format(self.submit_cmd) + temp_job_f_name + '\n'
            commands.append(comm)
            self.sims_done.add(new_sim_i)
        # process all jobs to run
        comm = '\ncd ' + sim_fd + '\n'
        commands.append(comm)
        for new_sim_i in sim.get('run'):
            new_sim = self.sim_set[new_sim_i]
            new_sim['flag_increase_count'] = False
            new_sim_info = self.process_sim(commands, new_sim, sim_fd_root, fnc2process, **kwargs)
            self.sims_done.add(new_sim_i)
        # write the job file with the commands
        #gen_slurm_job(job_f, commands, **kwargs) # job_f, commands, kwargs (e.g. slurm_kw)
        return job_f, sim_fd

    def generate_sim_files_jobs(self, sim_fd_root, fnc2process, **kwargs):
        submit_commands = []
        commands = []
        self.sims_done = set()
        for i, sim in enumerate(self.sim_set):
            if i not in self.sims_done:
                sim_info = self.process_sim(commands, sim, sim_fd_root, fnc2process, **kwargs)
                job_f, job_fd = sim_info
                self.sims_done.add(i)
                temp_job_f_name = os.path.split(job_f)[1]
                slurm_job_kwargs = sim.get('slurm_job_kwargs', kwargs)
                gen_slurm_job(job_f, commands, job_fd=job_fd, **slurm_job_kwargs) # job_f, commands, kwargs (e.g. slurm_kw)
                commands = []
                comm = 'cd ' + sim_info[1] + '\n{:} '.format(self.submit_cmd) + temp_job_f_name + '\n\n'
                print(comm)
                submit_commands.append(comm)
        return commands

SimulationSet._add_defaults(dict(submit_cmd='sbatch'), flag_set=1)

class Initial_FE(SimulationSet):
    """
    inital set of simulations for FE calculations
        generally consist of a set of shorter pre-equlibration simulations
        and a set of longer production FE calculations
        prepared for a set of lambda points
    """

    def __init__(self, t_eq, t_sim, LPs=None, flag_run_eq=True, max_t_sim=None, **kwargs):
        self.generate_sim_set(t_eq, t_sim, LPs=LPs, flag_run_eq=flag_run_eq, max_t_sim=max_t_sim, **kwargs)

    def generate_sim_set(self, t_eq, t_sim, LPs = None, flag_run_eq=True, max_t_sim=None, **kwargs):
        """
        gets a simulation set for FE calculations (dict) - saves it self.sim_set
        :param t_eq: equilibation time [ns]
        :param t_sim: production run time [ns]
        :param LPs: list of lambda points
        :param flag_run_eq: run all eq simulations in one job
        :param max_t_sim: max time for a signle simulation [ns] (e.g. time limit for slurm)
        :param kwargs:
            N_lam: number of equidistant lambda points (for LPs) - 11 by default
            flag_eq_sub: if defined as True, split eq simulations in different job files and submit separately
        :return:
            set of simulations (dictionary) - also saves it in self.sim_set
        """
        if LPs is None:
            N_lam = kwargs.get('N_lam', 11)
            LPs = get_lset(np.linspace(0, 1., N_lam))
        sim_set = []
        c = 0
        current_eq_sim_id = 0
        t_eq_total = 0
        if max_t_sim:
            max_t_sim = float(max_t_sim)
            eq_sim_Ts = self.split_single_sim(t_eq, max_t_sim)
            t_eq_total = eq_sim_Ts[0]
            FE_sim_Ts = self.split_single_sim(t_sim, max_t_sim)
        else:
            eq_sim_Ts = [t_eq]
            FE_sim_Ts = [t_sim]
        flag_eq_sub = kwargs.get('flag_eq_sub')
        for i,l in enumerate(LPs):
            if flag_run_eq:
                next_eq_sims = list(range(c+1, c+len(eq_sim_Ts)))
                if i != len(LPs) - 1:
                    next_eq_sims.append(c+1+len(FE_sim_Ts))
                else:
                    next_eq_sims.append(None)
                for eq_i, temp_t in enumerate(eq_sim_Ts):
                    subs = []
                    runs = []
                    # next eq simulation
                    next_eq_sim_i = next_eq_sims[eq_i]
                    if next_eq_sim_i is not None:
                        if flag_eq_sub:
                            subs.append(next_eq_sim_i)
                        else:
                            if max_t_sim:
                                next_eq_sim_t = eq_sim_Ts[(eq_i + 1) % len(eq_sim_Ts)]
                                t_eq_total += next_eq_sim_t
                                if t_eq_total > max_t_sim:
                                    subs.append(next_eq_sim_i)
                                    t_eq_total = next_eq_sim_t
                                else:
                                    runs.append(next_eq_sim_i)
                            else:
                                runs.append(next_eq_sim_i)
                    # next FE simulation
                    next_FE_sim = c+len(eq_sim_Ts)
                    subs.append(next_FE_sim)
                    new_sim = dict(eq=True, l=l, t=temp_t, sub=tuple(subs), run=tuple(runs))
                    sim_set.append(new_sim)
                    c += 1
            for FE_i, temp_t in enumerate(FE_sim_Ts):
                if FE_i != len(FE_sim_Ts) - 1:
                    subs = (c+1, )
                else:
                    subs = ()
                new_sim = dict(eq=False, l=l, t=temp_t, sub=subs, run=(), sys_name=self.FE_run_name)
                sim_set.append(new_sim)
                c+=1
        self.sim_set = sim_set

    def __format_lp(self, lp):
        return ('{:.'+str(self.LP_digits)+'f}').format(lp)

    def __get_sim_file_fd(self, sim, sim_fd_root, lp, flag_eq, **kwargs): # many of these def values should be encoded using the Defaults class!
        """
        based on the input simulation, generates a folder; job file name; and simulation name
        :param sim: simulation setup (usually a dict with setup params)
        :param sim_fd_root: folder in which all the new subfolders and job files will be created (creates it if doesn't exist)
        :param lp: lambda point
        :param flag_eq: flag that specifies if the simulation is a production or eq
        :param **kwargs:
            sims_count: count the simulations - dict with keys based on lambda points
            eq_lp_flag: count eq simulations based on lambda points individually or just normally (all count the same) - False by default
        :return:
            job_f: job file name
            sim_fd: folder in which this simulation will be done
            name: simulation name (this will also define the in/out file names)
        """
        if 'sims_count' in kwargs:
            self.sims_count = kwargs['sims_count']
        if not hasattr(self, 'sims_count'):
            self.sims_count = {'eq':{}, 'FE':{}}
        sims_count = self.sims_count
        lp = fix_float(lp)
        eq_lp_flag = kwargs.get('eq_lp_flag')
        if flag_eq:
            sim_key = 'eq'
            sim_fd = os.path.join(sim_fd_root, self.eq_fd)
            if eq_lp_flag:
                sim_key = lp
            sim_count = self.get_sim_count(sims_count['eq'], sim_key, 
                                            flag_increase_count=sim.get('flag_increase_count', True))
            job_f = os.path.join(sim_fd,  '{:}_{:}.sh'.format(self.eq_run_file, sim_count))
            name = '{:}_{:}'.format(self.eq_name, self.__format_lp(lp))
        else:
            sim_fd = os.path.join(sim_fd_root, '{:}_{:}'.format(self.FE_run_fd, self.__format_lp(lp)))
            sim_count = self.get_sim_count(sims_count['FE'], lp)
            job_f = os.path.join(sim_fd, '{:}_{:}.sh'.format(self.FE_run_file, sim_count))
            name = '{:}_{:}'.format(self.FE_run_name, sim_count)
        if not os.path.isdir(sim_fd):
            os.mkdir(sim_fd)
        return job_f, sim_fd, name

    def get_sim_file_fd(self, sim, sim_fd_root, **kwargs):
        return self.__get_sim_file_fd(sim, sim_fd_root, sim['l'], sim['eq'], **kwargs)

_Initial_FE_defs = dict(eq_fd='pre_eq', eq_run_file='run_pre_eq', eq_name='preeq_l',
                        FE_run_fd='L', FE_run_file='run_FE', FE_run_name='FE_MD', LP_digits=3)
Initial_FE._add_defaults(_Initial_FE_defs, flag_set=1)

# for compatibility with the old scripts...
def run_intial_FE(t_eq, t_sim, LPs=None, flag_run_eq=True, max_t_sim=None, **kwargs):
    """
    gets a simulation set for FE calculations
    :param t_eq: equilibation time [ns]
    :param t_sim: production run time [ns]
    :param LPs: list of lambda points
    :param flag_run_eq: run all eq simulations in one job
    :param max_t_sim: max time for a signle simulation [ns] (e.g. time limit for slurm)
    :param kwargs:
        N_lam: number of equidistant lambda points (for LPs) - 11 by default
        flag_eq_sub: if defined as True, split eq simulations in different job files and submit separately
    :return:
        set of simulations (dictionary)
    """
    sim_set = Initial_FE(t_eq, t_sim, LPs=LPs, flag_run_eq=flag_run_eq, max_t_sim=max_t_sim, **kwargs)
    return sim_set.sim_set

