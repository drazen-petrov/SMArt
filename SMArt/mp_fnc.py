from SMArt.incl import os, subprocess, Pool, do_warn

def run_n_mp(n, f2run, args):
    p = Pool(n)
    res = None
    try:
        res = p.map_async(f2run, args).get()
    except:
        print("Unexpected error in f2run:\n")
    p.terminate()
    return res


def _change_fd(comm_fd):
    if isinstance(comm_fd, str):
        return comm_fd
    elif len(comm_fd)==2 and isinstance(comm_fd[0], str) and isinstance(comm_fd[1], str):
        if os.path.isdir(comm_fd[1]):
            os.chdir(comm_fd[1])
            return comm_fd[0]
        else:
            do_warn('folder not found ' + comm_fd[1])
    else:
        do_warn('wrong input format\n(command, fd) expected\ngot: ' + str(comm_fd))

def _run_comm_v1(comm_fd):
    comm = _change_fd(comm_fd)
    if comm:
        os.system(comm)

def _run_comm_v2(comm_fd, flag_print=False, **kwargs):
    """
    :param comm_fd: either str (just command), or [str, str] command, folder
    :param flag_print:
    :param kwargs:
        shell = True by default
        timeout
        other from subprocess.check_output
    :return:
    """
    kwargs = dict(kwargs)
    kwargs['shell'] = kwargs.get('shell', True)
    comm = _change_fd(comm_fd)
    if comm:
        try:
            s = subprocess.check_output(comm, **kwargs)
            if flag_print:
                print(s.decode("utf-8"))
        except subprocess.TimeoutExpired:
            do_warn('timed out after {:}'.format(kwargs.get('timeout')))

def run_comm_n_mp(n, comm_fd, fnc_v = 2):
    run_comm_fnc = {1:_run_comm_v1, 2:_run_comm_v2}[fnc_v]
    run_n_mp(n, run_comm_fnc, comm_fd)
