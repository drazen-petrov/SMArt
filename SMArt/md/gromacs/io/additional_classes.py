from .gm_io_defaults import gm_io_Defaults

class Define(gm_io_Defaults):
    """gromacs define class
    :var
        id
        v - values
        s - string (line)"""
    container2write = 'defines'

    def __init__(self, split_line = None, line = None, flag_generate_line = True, **kwargs):
        if split_line is None and line is None:
            raise Exception('split_line is None and line is None')
        if split_line is None:
            split_line = line.split()
        if self._define_directive == split_line[0]:
            split_line = split_line[1:]
        self.id = split_line[0]
        self.v = split_line[1:]
        if line is None and flag_generate_line:
            line = self.__generate_line()
        self.s = line

    def __str__(self):
        return str(self.v)

    def __repr__(self):
        return str(self.v)

    def __generate_line(self):
        txt = self._define_directive + ' {:}' + ' {:>9}' * len(self.v)
        temp = [str(self.id)]
        for i in self.v:
            temp.append(str(i))
        return txt.format(*temp) + '\n'

    def write_define(self, from_str = False, **kwargs):
        if from_str:
            if hasattr(self, 's'):
                return self.s
        return self.__generate_line()


class ParsedSegments(gm_io_Defaults):
    """
    each parsed segment is an instance of this class
    ensures that writing gets done properly
    :attr
        seg_type (include or define)
        seg_name
        curr_dir - current directive (e.g. [ bonds ])
        curr_inc - current include (included file)
        curr_if - current if statement
        curr_seg_num
        set_curr_stack
    """

    container2write = '_segments'

    def __init__(self, seg_kw = None, seg_name = None, curr_dir = None, curr_inc = None, curr_if = None,
                 curr_seg_num = None,  **kwargs):
        """
        :param seg_kw: defines seg_type, e.g. #include or #define
        :param seg_name:
        :param curr_dir:
        :param curr_inc:
        :param curr_if:
        :param curr_seg_num:
        :param kwargs:
        """
        self.seg_type = self._get_segtype(seg_kw, kwargs.get('seg_type')) # include, define, undef, gm_directive, if
        self.seg_name = seg_name
        self.curr_dir = curr_dir
        self.curr_inc = curr_inc
        self.curr_if = curr_if
        self.s = kwargs.get('s', str())
        self.curr_seg_num = curr_seg_num
        curr_stack = kwargs.get('curr_stack', None)
        if curr_stack is not None:
            self.set_curr_stack(curr_stack)
            curr_stack[3]+=1
        file_fd_path = kwargs.get('file_fd_path')
        if file_fd_path:
            self.fd, self.f_name, self.f_path = file_fd_path

    def __str__(self):
        return str((self.curr_seg_num, self.seg_type, self.seg_name))

    def __repr__(self):
        return self.__str__()

    def set_curr_stack(self, curr_stack):
        # sets the stack, self.curr_dir, self.curr_inc, self.curr_if, self.curr_seg_num
        self.curr_dir = curr_stack[0]
        self.curr_inc = tuple(curr_stack[1]) # this is needed to avoid changes
        self.curr_if = tuple(curr_stack[2]) # this is needed to avoid changes
        self.curr_seg_num = curr_stack[3]

    def get_curr_stack(self):
        # gets the stack, self.curr_dir, self.curr_inc, self.curr_if, self.curr_seg_num
        return (self.curr_dir, self.curr_inc, self.curr_if, self.curr_seg_num)

    @property
    def gm_ord(self):
        return self.get_curr_stack()

    def _write_include(self, abspath = False):
        if abspath:
            return self._include_directive + ' ' + self.f_path
        else:
            return self.s

    def _write_define(self, **kwargs): # can be different from write_define from Define!
        return self.s

    def _write_undef(self, **kwargs):
        return self.s

    def _write_if(self, **kwargs):
        return self.s

