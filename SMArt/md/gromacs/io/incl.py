from SMArt.incl import os, copy, OrderedDict, do_warn
from SMArt.incl import Defaults, GeneralContainer, FileStream, StringStream, split_gen, get_id_string
from SMArt.md.incl import *
#from SMArt.md.incl import _gr_title_gm_system, _file_type_ext_map, _add_DescriptionPart
from .gm_io_defaults import gm_io_Defaults, gm_FDs_Defaults, gm_FF_Defaults
from .additional_classes import Define, ParsedSegments


class GromacsStream(gm_FDs_Defaults, gm_io_Defaults):
    def __just_for_readability(self):
        self._gm_path_FDs = []
        self.f = None
        self._gm_eol = None
        self._gm_comm_indicators = None

    @property
    def _comm_indicators(self):
        return self._get_gm_comm_indicators()

    def lines(self):
        """this is implemented in GromacsFile or GromacsString"""
        yield None

    def _get_fd(self):
        """this is implemented in GromacsFile or GromacsString"""
        pass

    @staticmethod
    def _rem_q(f_name):  # remove quotes
        qs = ["\"", "\'"]
        if f_name[0] in qs:
            f_name = f_name[1:]
        if f_name[-1] in qs:
            f_name = f_name[:-1]
        return f_name

    def _find_file(self, f_path, fd=False):
        if not isinstance(f_path, str):
            f_path = str(f_path)
        if fd:
            temp_path = os.path.join(fd, f_path)
            if os.path.isfile(temp_path):
                return temp_path
            else:
                raise Exception("file doesn't exist: " + temp_path)
        if os.path.isfile(f_path) and os.path.abspath(f_path)==f_path:
            return f_path
        temp_path = os.path.join(self._get_fd(), f_path)
        if os.path.isfile(temp_path):
            return temp_path
        else:
            temp_path = self._find_file_gm_path(f_path)
            if temp_path:
                return temp_path
            raise Exception("file doesn't exist: " + f_path)

    def _open_file(self, f_path, fd=False):
        # searches and opens a file for the file in the local folder and in FF_DB folder - returns error if not found
        f_path = self._find_file(f_path, fd)
        fd, file_name = os.path.split(f_path)
        return GromacsFile(f_path), fd, file_name

    @staticmethod
    def _check_flag_w(if_stack):
        """makes an if_stack - makes sure to read nested if statements"""
        for i in if_stack:
            if not i:
                return 0
        return 1

    @staticmethod
    def read_define(temp, defines, line = None, temp_seg = None):  # reads difine lines - interaction types (e.g. gd_1)
        temp_def = Define(temp, line)
        #temp_def.gm_ord = temp_seg
        defines[temp[1]] = temp_def

    @staticmethod
    def read_undef(temp, defines):
        if temp[1] in defines:
            del(defines[temp[1]])

    @staticmethod
    def find_directive(l):
        if "[" in l and "]" in l:
            return l[l.find("[") + 1:l.find("]")].strip()
        else:
            return False

    def _read_lines_generator(self, defines = None, curr_stack = None):
        "generator that return line by line - connects lines that have to be connected and takes care of if statements"
        if curr_stack is None:
            curr_stack = [None, list(), list(), 0]
        if defines is None:
            defines = OrderedDict()
        flag_w = 1
        line = ''
        if_stack = []
        for temp_l in self.lines():
            if temp_l:
                line += temp_l
            if line.strip() == '':
                line = ''
                continue
            line, flag_eol = self._remove_gm_eol(line)
            if flag_eol:
                continue
            if line and not line.endswith('\n'):
                line+='\n'
            temp = line.split()
            for if_seg in curr_stack[2]:
                if_seg.s += line
            if temp[0] in self._gm_kw:
                if temp[0] == self._include_directive:
                    if flag_w:
                        new_f_name = self._rem_q(temp[1])  # because it comes with " "
                        gms, fd, file_name = self._open_file(new_f_name)
                        temp_seg = ParsedSegments(temp[0], curr_stack = curr_stack, s = line,
                                                  file_fd_path = (fd, file_name, gms.f_path))
                        #temp_seg.fd, temp_seg.f_name, temp_seg.f_path = fd, file_name, gms.f_path
                        curr_stack[1].append(temp_seg)
                        #curr_stack[3] += 1
                        yield temp_seg
                        for included_line in gms._read_lines_generator(defines, curr_stack = curr_stack):
                            yield included_line
                        gms._remove_f()
                        _ = curr_stack[1].pop()
                        curr_stack[3] += 1
                elif temp[0] == "#define":
                    if flag_w:
                        temp_seg = ParsedSegments(temp[0], seg_name=temp[1], curr_stack=curr_stack)
                        temp_seg.s = line
                        #curr_stack[3] += 1
                        self.read_define(temp, defines, line, temp_seg)
                        yield temp_seg
                elif temp[0] == "#undef":
                    if flag_w:
                        temp_seg = ParsedSegments(temp[0], curr_stack=curr_stack)
                        temp_seg.s = line
                        #curr_stack[3] += 1
                        yield temp_seg
                        self.read_undef(temp, defines)
                elif temp[0] == "#ifdef":
                    temp_seg = ParsedSegments(temp[0], curr_stack = curr_stack)
                    temp_seg.s = line
                    #curr_stack[3] += 1
                    curr_stack[2].append(temp_seg)
                    yield temp_seg
                    if temp[1] in defines:
                        if_stack.append(True)
                    else:
                        if_stack.append(False)
                        flag_w = 0
                elif temp[0] == "#ifndef":
                    temp_seg = ParsedSegments(temp[0], curr_stack=curr_stack)
                    temp_seg.s = line
                    #curr_stack[3] += 1
                    curr_stack[2].append(temp_seg)
                    yield temp_seg
                    if not temp[1] in defines:
                        if_stack.append(True)
                    else:
                        if_stack.append(False)
                        flag_w = 0
                elif temp[0] == "#else":
                    if_stack[-1] = not if_stack[-1]
                    flag_w = self._check_flag_w(if_stack)
                elif temp[0] == "#endif":
                    _ = if_stack.pop()
                    curr_stack[3] += 1
                    _ = curr_stack[2].pop()
                    flag_w = self._check_flag_w(if_stack)
                else:
                    raise Exception
            elif flag_w:
                temp_directive = self.find_directive(temp_l)
                if temp_directive:
                    temp_seg = ParsedSegments(seg_type='gm_directive', seg_name=temp_directive, curr_stack = curr_stack)
                    temp_seg.s = line
                    #curr_stack[3] += 1
                    curr_stack[0] = temp_seg
                    yield temp_seg
                else:
                    yield line
            line = ""

    def read_lines(self, defines = None, curr_stack = None):
        if curr_stack is None:
            curr_stack = [None, list(), list(), 0]
        if defines is None:
            defines = OrderedDict()
#        if hasattr(self, 'defines'):
#            self.defines.update(defines)
#        else:
#            self.defines = defines
        self.defines = defines
        self._buffer_l = self._read_lines_generator(self.defines, curr_stack)
        for i in self._buffer_l:
            yield i
        self._buffer_l = ''

    def next_line(self):
        if self._buffer_l:
            return next(self._buffer_l)
        return ''

    def write(self, s2w):
        self.f.write(s2w)
        if isinstance(self, GromacsString):
            self.s+=s2w

    def write_lines(self, lines):
        s2w = ''
        for l in lines:
            if not l.endswith('\n'):
                l+='\n'
            s2w+=l
        self.write(s2w)


class GromacsFile(FileStream, GromacsStream):
    """handler class for gromos files - reading / writing functions"""
    def _get_fd(self):
        return os.path.split(self.f_path)[0]

class GromacsString(StringStream, GromacsStream):
    """handler class for gromos string - reading / writing functions"""
    def _get_fd(self):
        return os.getcwd()


direction_interaction_map = OrderedDict()
direction_interaction_map['bonds'] = BondType
direction_interaction_map['pairs'] = PairType
direction_interaction_map['pairs_nb'] = PairNBType
direction_interaction_map['angles'] = AngleType
direction_interaction_map['dihedrals'] = DihedralType
direction_interaction_map['impropers'] = ImproperType
# this is not for reading - only for writing (see below reverse_interaction_map)
direction_interaction_map['cmap'] = cmap
direction_interaction_map['constraints'] = ConstraintType
direction_interaction_map['settles'] = SettleType
direction_interaction_map['virtual_sites2'] = VirtualSite2Type
direction_interaction_map['virtual_sites3'] = VirtualSite3Type
direction_interaction_map['virtual_sites4'] = VirtualSite4Type
direction_interaction_map['virtual_sitesn'] = VirtualSitenType
direction_interaction_map['position_restraints'] = Position_r
direction_interaction_map['distance_restraints'] = Distance_r
direction_interaction_map['orientation_restraints'] = Orientation_r
direction_interaction_map['angle_restraints'] = Angle_r
direction_interaction_map['angle_restraints_z'] = Angle_r_z

direction_interaction_type_map = OrderedDict()
direction_interaction_type_map['nonbond_params'] = vdWType
direction_interaction_type_map['pairtypes'] = PairType
direction_interaction_type_map['bondtypes'] = BondType
direction_interaction_type_map['angletypes'] = AngleType
direction_interaction_type_map['dihedraltypes'] = DihedralType
direction_interaction_type_map['impropertypes'] = ImproperType
direction_interaction_type_map['cmaptypes'] = cmap
for temp_dir in direction_interaction_map:
    if temp_dir.endswith('s') and temp_dir[:-1] + 'types' not in direction_interaction_type_map:
        direction_interaction_type_map[temp_dir[:-1] + 'types'] = direction_interaction_map[temp_dir]

directive_name_map = {'impropers':'dihedrals', 'impropertypes':'dihedraltypes'}

ff_directives = [['defaults'], ['atomtypes'], ['nonbond_params', 'pairtypes']]
ff_directives.append([])
for d in direction_interaction_type_map:
    if d not in ff_directives[-2]:
        ff_directives[-1].append(d)

#top_directives = copy.deepcopy(ff_directives)
#top_directives.extend([['moleculetype'], ['system'], ['molecules']])
top_directives = [['system'], ['molecules']]


class GromacsParser(GeneralContainer, Defaults):
    """Genearal parsing functions

    interesting attrib
        gm_curr_stack - defines stack of directives, defines, IFs, etc.
                        self.curr_dir, self.curr_inc, self.curr_if, self.curr_seg_num
        _gm_curr_fnc
        _gm_curr_int_type
        undefined_directives"""

    def __just_for_readability(self):
        self.gm_ff_defaults = list()
        self.__parse_gm_v = None

    Define = Define

    def add_define(self, split_line = None, line = None, flag_generate_line = True, **kwargs):
        """
        :param kwargs:
            split_line = None, line = None, flag_generate_line = True
            add2container kwargs
        :return:
        """
        kwargs['split_line'] = split_line
        kwargs['line'] = line
        kwargs['flag_generate_line'] = flag_generate_line
        kwargs['create'] = True
        self.add2container(Define(**kwargs), **kwargs)

    @staticmethod
    def __parse_from_source(parsing_from, parse_from_file=True):
        if isinstance(parsing_from, (GromacsFile, GromacsString)):
            return parsing_from
        if parse_from_file:
            parsing_from = GromacsFile._find_file_gm_path_cl(parsing_from)
            return GromacsFile(parsing_from)
        return GromacsString(parsing_from)

    @staticmethod
    def exchange_defs(l_split, defs):
        temp = []
        for i in l_split:
            if i in defs:
                temp.extend(defs[i])
            else:
                temp.append(i)
        return temp

    def _gm_find_parse_fnc(self, directive_name, fnc_pref = None, fnc_suf = None, **kwargs):
        if fnc_pref is None:
            fnc_pref=''
        if fnc_suf is None:
            fnc_suf = ''
        fnc_directive_name = fnc_pref + directive_name + fnc_suf
        if fnc_directive_name!=directive_name:
            if hasattr(self, '_' + fnc_directive_name + '_parser'):
                fnc_name = getattr(self, '_' + fnc_directive_name + '_parser')
                return getattr(self, fnc_name), False
        if hasattr(self, '_' + directive_name + '_parser'):
            fnc_name = getattr(self, '_' + directive_name + '_parser')
            return getattr(self, fnc_name), False
        if not kwargs.get('ignore_active_objects'):
            if hasattr(self, '_curr_bbmol'):
                temp_f = self._curr_bbmol._gm_find_parse_fnc(directive_name, fnc_pref = fnc_pref, fnc_suf = fnc_suf,
                                                             allow_gm_find_parse_fnc_not_found = 1, **kwargs)
                if temp_f:
                    if directive_name in direction_interaction_map:
                        self._curr_bbmol._gm_curr_int_type = direction_interaction_map[directive_name]
                    return temp_f[0], self._curr_bbmol
        if kwargs.get('allow_gm_find_parse_fnc_not_found'):
            return
        else:
            do_warn('unknown directive read in as text: ' + directive_name)
            if directive_name not in self.undefined_directives:
                self.undefined_directives[directive_name] = []
            self.undefined_directives[directive_name].append(self.gm_curr_stack[0])
            return self.__read_unknown_directive, False

    def __read_unknown_directive(self, l, *args, **kwargs):
        self.gm_curr_stack[0].s += l

    def _add_gm_ord(self, item, **kwargs):
        # sets gm_ord - curr_dir, curr_inc, curr_if, curr_seg_num
        temp_gm_curr_stack = kwargs.get('gm_curr_stack', None)
        if temp_gm_curr_stack is None:
            temp_gm_curr_stack = self.gm_curr_stack
        item.gm_ord = (temp_gm_curr_stack[0], tuple(temp_gm_curr_stack[1]),
                       tuple(temp_gm_curr_stack[2]), temp_gm_curr_stack[3])


    def __parse_gm_v1(self, parsing_from, parse_from_file, **kwargs):
        if not hasattr(self, 'gm_curr_stack'):
            curr_stack = kwargs.get('gm_curr_stack', [None, list(), list(), 0])
            self.gm_curr_stack = curr_stack
            temp_seg = ParsedSegments(curr_stack = self.gm_curr_stack)
            curr_stack[0] = temp_seg
            #curr_stack[3] += 1
            self.add2container(temp_seg, create=True, db_type=list)
        else:
            curr_stack = self.gm_curr_stack
        kwargs['gm_curr_stack'] = curr_stack
        if not hasattr(self, 'defines'):
            self.defines = OrderedDict()
        defines = kwargs.get('defines')
        if defines is not None:
            self.defines.update(defines)
        defines = self.defines
        gms = self.__parse_from_source(parsing_from, parse_from_file)
        temp_title = DescriptionPart(form = 'gm', **kwargs)
        temp_title._add_source(gms.f_path, additional_txt = 'GROMACS')
        self._add_gm_ord(temp_title)
        self.add2container(temp_title, create=True, db_type=list)
        if not hasattr(self, 'parsed_from'):
            self.parsed_from = []
        self.parsed_from.append(gms)
        if not hasattr(self, 'undefined_directives'):
            self.undefined_directives = OrderedDict()
        self._gm_curr_fnc = self.__read_unknown_directive
        for l in gms.read_lines(defines = defines, curr_stack = self.gm_curr_stack):
            if not isinstance(l, ParsedSegments):
                self._gm_curr_fnc(l, defines, **kwargs)
            else:
                self.add2container(l)
                if l.seg_type=='gm_directive':
                    if l.seg_name in direction_interaction_map:
                        self._gm_curr_int_type = direction_interaction_map[l.seg_name]
                    if l.seg_name in direction_interaction_type_map:
                        self._gm_curr_ff_int_type = direction_interaction_type_map[l.seg_name]
                    self._gm_curr_fnc, parsing_object = self._gm_find_parse_fnc(l.seg_name, **kwargs)
                elif l.seg_type=='gm_define':
                    self._add_gm_ord(self.defines[l.seg_name])
                elif l.seg_type == 'include':
                    temp_title = DescriptionPart(form='gm', **kwargs)
                    temp_title._add_source(l.f_path, additional_txt='GROMACS')
                    self._add_gm_ord(temp_title)
                    self.add2container(temp_title, create=True, db_type=list)
        if not kwargs.get('ignore_title_lines'):
            seg0_lines = self._segments[0].s.split('\n')
            if seg0_lines:
                self.get_container(DescriptionPart.container2write)[0].lines.extend(seg0_lines)
        del(self._gm_curr_fnc)
        gms._remove_f()

    def parse_gm(self, parse_from, parse_from_file = True, **kwargs):
        """
        parsing gromacs format
        :param parse_from:
        :param parse_from_file:
        :param kwargs:
            gm_curr_stack
            defines - gromacs defines
            ignore_title_lines
        :return:
        """
#        temp_fd = os.getcwd()
        temp_f = getattr(self, self.__parse_gm_v)
        temp_f(parse_from, parse_from_file, **kwargs)
#        os.chdir(temp_fd)

GromacsParser._add_defaults({'__parse_gm_v':'__parse_gm_v1'}, flag_set=True)


class GromacsWriter(GeneralContainer, Defaults):
    """Genearal writing functions"""
    def __just_for_readability(self):
        self.defines = OrderedDict()

    def __write_directive_line_v1(self, gs, directive, **kwargs):
        if directive in directive_name_map:
            directive = directive_name_map[directive]
        gs.write('\n[ {:} ]\n'.format(directive))

    def _write_directive_line(self, gs, directive, **kwargs):
        temp_f = getattr(self, self.__write_directive_line_v)
        temp_f(gs, directive, **kwargs)

    def __write_define_v1(self, gs, *defs, **kwargs):
        temp_defines = kwargs.get('defines2write', {})
        if temp_defines:
            if kwargs.get('use_self_defines', True):
                self_defs = getattr(self, 'defines', {})
                if self_defs:
                    temp_defines = copy.deepcopy(temp_defines)
                    temp_defines.update(self_defs)
        else:
            temp_defines = getattr(self, 'defines', {})
        if defs:
            for d in defs:
                if self._check_flag_write(temp_defines[d], **kwargs):
                    gs.write(temp_defines[d].write_define(**kwargs))
        else:
            gm_defs_from_seg = set()
            if kwargs.get('flag_segment_order'):
                if not kwargs.get('flag_define_noseg', True):
                    return
                gm_parsed_segments = kwargs.get('gm_parsed_segments', getattr(self, '_segments', list()))
                for seg in gm_parsed_segments:
                    if seg.seg_type == 'gm_define':
                        gm_defs_from_seg.add(seg.seg_name)
            for d in temp_defines:
                if d not in gm_defs_from_seg:
                    if kwargs.get('flag_segment_order') and hasattr(temp_defines[d],'gm_ord'):
                        continue
                    if self._check_flag_write(temp_defines[d], **kwargs):
                        gs.write(temp_defines[d].write_define(**kwargs))

    def _write_define(self, gs, **kwargs):
        """
        :param gs:
        :param kwargs:
            defs - list of defines to write
        :return:
        """
        defs = kwargs.get('defs', [])
        temp_f = getattr(self, self.__write_define_v)
        temp_f(gs, *defs, **kwargs)

    def __write_TITLE_v1(self, gs, **kwargs):
        flag_title = kwargs.get('flag_title')
        temp_lines = []
        first_line = []
        temp_cont = self.get_container(DescriptionPart, flag_class=True, allow_not_found = 1)
        if temp_cont is None:
            return
        if flag_title != -1:
            temp_cont = temp_cont[:flag_title]
        for temp_title in temp_cont:
            temp_lines.extend(first_line)
            first_line = ('', '')
            for temp_l in temp_title.lines:
                temp_lines.append(gs._comm_indicators[0] + temp_l)
        temp_lines.append('')
        gs.write_lines(temp_lines)

    def _write_gm_TITLE(self, gs, **kwargs):
        temp_f = getattr(self, self.__write_TITLE_v)
        temp_f(gs, **kwargs)

    def _write_include(self, gs, include_segment_path, **kwargs):
        """
        :param gs:
        :param include_segment_path:
        :param kwargs:
            full_path
            file_name
            from_str
        :return:
        """
        path2write = None
        if isinstance(include_segment_path, str):
            path2write = include_segment_path
        elif kwargs.get('full_path'):
            path2write = include_segment_path.f_path
        elif kwargs.get('file_name'):
            path2write = include_segment_path.f_name
        elif kwargs.get('from_str'):
            gs.write(include_segment_path.s)
        else:
            path2write = include_segment_path.f_path
        if path2write:
            write_format = '{0:} "{1:}"\n'
            gs.write(write_format.format(gs._gm_kw[0], path2write))
            #gs.write(gs._gm_kw[0] + ' "' + path2write + '"\n')

    def _write_if(self, gs, temp_if, **kwargs):
        gs.write(temp_if.s)

    def write_gromacs_format(self, gs, directives = None, flag_include = False, flag_if = False,
                                  flag_defines = True, flag_title = 1, **kwargs):
        directives = self.__copy_directives(directives)
        temp_f = getattr(self, self.__write_gromacs_format_v)
        full_kwargs = self.__update_kwargs4write_gromacs_format(flag_include = flag_include, flag_if = flag_if,
                                  flag_defines = flag_defines, flag_title = flag_title, **kwargs)
        temp_f(gs, directives = directives, **full_kwargs)
        if kwargs.get('flag_close'):
            gs.f.close()

    @staticmethod
    def __copy_directives(directives):
        if directives is None:
            directives = []
        flag = 1
        for grp_dir in directives:
            if not hasattr(grp_dir, '__iter__') or isinstance(grp_dir, str):
                flag = 0
                break
        new_directives = []
        if flag:
            for grp_dir in directives:
                new_directives.append(list(grp_dir))
        else:
            new_directives.append(list(directives))
        return new_directives

    @staticmethod
    def __update_kwargs4write_gromacs_format(**kwargs):
        default_KWs = dict(flag_title=True, flag_include=False, flag_if=False, flag_defines=True,
                           flag_use_define = True, flag_ff_defines = True)
        for temp_kw in default_KWs:
            kwargs[temp_kw] = kwargs.get(temp_kw, default_KWs[temp_kw])
        return kwargs

    def __set_gm_writing_state(self):
        if not hasattr(self, '_GromacsWriter__writing_state_var'):
            self.__writing_state_var = [([], [], [], [], [])] # gs, include, if, defines, directives

    @property
    def _gm_writing_state(self):
        self.__set_gm_writing_state()
        if not self.__writing_state_var:
            self._add_gs2writing_state(GromacsString())
        return self.__writing_state_var

    def _add_gs2writing_state(self, gs):
        self.__set_gm_writing_state()
        self.__writing_state_var[-1][0].append(gs)

    def _add_include_gs2writing_state(self, f_path, **kwargs):
        gs = self._gm_writing_state[-1][0][-1]
        self._write_include(gs, f_path, **kwargs)
        gs = GromacsFile(f_path, True)
        self._add_gs2writing_state(gs)
        return gs

    def _remove_gs_from_writing_state(self):
        old_gs = self._gm_writing_state[-1][0].pop(-1)
        old_gs.f.close()
        return self._gm_writing_state[-1][0][-1]

    def _add_state2writing_state(self, gs):
        self.__set_gm_writing_state()
        self.__writing_state_var.append(([gs], [], [], [], [])) # gs, include, if, defines, directives

    @staticmethod
    def __transform_directives_sep2include_files(directives, **kwargs):
        """
        :param directives:
        :param sep2include_files:
        :param kwargs:
            last_inc_file = None
        :return:
        """
        sep2include_files = kwargs.get('sep2include_files', [])
        if not sep2include_files:
            sep2include_files = [kwargs.get('last_inc_file', None)]
        flag_transform = 0
        for i in sep2include_files:
            if not isinstance(i, tuple) or len(i)!=2:
                flag_transform = 1
                break
            else:
                assert (i[0] is None or isinstance(i[0], str)) and isinstance(i[1], int) and i[1]>0
        if flag_transform:
            sep2include_files_transformed = []
            c_curr = 1
            curr_directive = sep2include_files.pop(0)
            while sep2include_files:
                next_directive = sep2include_files.pop(0)
                if curr_directive == next_directive:
                    c_curr += 1
                else:
                    sep2include_files_transformed.append((curr_directive, c_curr))
                    c_curr = 1
                    curr_directive = next_directive
            sep2include_files_transformed.append((curr_directive, c_curr))
        c = 0
        temp_dir = []
        directives_transformed = []
        for sep_inc_file in sep2include_files:
            if sep_inc_file[0] is not None:
                assert sep_inc_file[0] not in temp_dir, 'non-consecutive usage of the same file name'
                temp_dir.append(sep_inc_file[0])
            new_c = c + sep_inc_file[1]
            directives_transformed.append((sep_inc_file[0], directives[c:new_c]))
            c = new_c
        directives_transformed.append((kwargs.get('last_inc_file', None), directives[c:]))
        return directives_transformed

    def __write_gromacs_format_v1(self, gs, directives = None, **kwargs):
        """
        :param gs:
        :param directives: [[dir1_1, dir1_2, ...], [dir2_1, dir2_2, ...], ...] in groups
        :param kwargs:
            flag_include
            flag_if
            flag_defines
            flag_use_define
            flag_title
            flag_segment_order
            check_self_flag_write
            flag_ff_defines
            sep2include_files = [None, None, f_path1, f_path1, None, f_path2, ...]
                                [(None, 2), (f_path1,2), (None, 1), (f_path2, 1), ...]
            last_inc_file = None
            flag_defines_first_include = False
        """
        assert hasattr(directives, '__iter__')
        self.__writing_state_var = kwargs.get('gm_writing_state', [([gs], [], [], [], [])])
        assert gs == self.__writing_state_var[-1][0][-1], 'gromacs stream != from __writing_state_var defined gs'
        self._segments_writing_stack = kwargs.get('gm_segments_writing_stack', list(getattr(self, '_segments',[])))
        if kwargs['flag_include']:
            kwargs['flag_segment_order'] = kwargs.get('flag_segment_order', True)
        if kwargs.get('check_self_flag_write'):
            flag_write = self._check_flag_write(self, **kwargs)
            if not flag_write:
                del(self.__writing_state_var)
                return
        directives_transformed = self.__transform_directives_sep2include_files(directives, **kwargs)
        self._write_gm_TITLE(gs, **kwargs)
        if kwargs['flag_defines']:############################################################### fix to write gm defaults
            if not kwargs.get('flag_defines_first_include'):
                self._write_define(gs, **kwargs)
#        for temp_fpath, grp_dir in zip(sep2include_files, directives):
        for temp_fpath, directives in directives_transformed:
            if temp_fpath:
                gs = self._add_include_gs2writing_state(temp_fpath, **kwargs)
                if kwargs.get('flag_defines_first_include') and kwargs['flag_defines']:
                    self._write_define(gs, **kwargs)
                    kwargs['flag_defines_first_include'] = False
            for grp_dir in directives:
                if kwargs.get('flag_segment_order'):
                    self.__write_gm_order_directives(gs, grp_dir, **kwargs)
                for temp_dir in grp_dir:
                    if kwargs.get('flag_segment_order'):
                        kwargs['flag_check_directive_inWS'] = 1
                    else:
                        kwargs['flag_check_directive_inWS'] = 0
                    self._write_directive(gs, temp_dir, **kwargs)
            if temp_fpath:
                gs = self._remove_gs_from_writing_state()
#                self._gm_writing_state[-1][0].pop(-1)
#                gs = self._gm_writing_state[-1][0][-1]
        if kwargs.get('flag_del_writing_state', True):
            del(self.__writing_state_var)

    def __write_gm_order_directives(self, gs, gm_grp_directives, **kwargs):
        if not self._segments_writing_stack:
            return
        segs2remove = []
        segs2remove_map = {}
        for i, seg in enumerate(self._segments_writing_stack):
            if seg.seg_type == 'gm_directive' and seg.seg_name not in gm_grp_directives:
                break
            segs2remove.append(i)
            segs2remove_map[seg] = i
            self.__write_gm_order_segment(gs, seg, **kwargs)
        for temp_seg in seg.gm_ord[1]:
            if temp_seg in segs2remove_map:
                del(segs2remove_map[temp_seg])
        for temp_seg in seg.gm_ord[2]:
            if temp_seg in segs2remove_map:
                del(segs2remove_map[temp_seg])
        for i in sorted(segs2remove_map.values(),reverse=True):
            seg = self._segments_writing_stack.pop(i)
            if seg.seg_type in ('if', 'include'):
                flag_write = self._check_flag_write(seg, **kwargs)
                if not flag_write:continue
                if seg.seg_type == 'if' and kwargs.get('flag_if'):
                    self._check_flag_if(seg, flag_if_temp_gm_object=True, **kwargs)
                if seg.seg_type == 'include' and kwargs.get('flag_include'):
                    self._check_flag_include(seg, flag_include_temp_gm_object=True, **kwargs)

    def __write_gm_order_segment(self, gs, gm_segment, **kwargs):
        if gm_segment.seg_type == 'gm_define':
            if kwargs['flag_defines']:
                if self._check_flag_write(gm_segment, **kwargs):
                    if kwargs.get('def_from_segment'):
                        gs.write(gm_segment.s)
                    else:
                        self._write_define(gs, defs = [gm_segment.seg_name], **kwargs)
        elif gm_segment.seg_type == 'gm_directive':
            kwargs['curr_directive'] = gm_segment
            #self._write_directive(gs, gm_segment.seg_name, curr_directive = gm_segment, **kwargs)
            self._write_directive(gs, gm_segment.seg_name, **kwargs)
            self._gm_writing_state[-1][4].append(gm_segment)

    def __write_directive_v1(self, gs, temp_dir, **kwargs):
        if temp_dir in direction_interaction_type_map:
            int_type = direction_interaction_type_map[temp_dir]
            temp_cont = self.get_intDB(int_type, cont_pref=self._gm_cont_pref, create=False, allow_not_found=1)
            if not temp_cont:return
            self._gm_curr_write_ff_int_type = int_type
        if temp_dir in direction_interaction_map:
            int_type = direction_interaction_map[temp_dir]
            self._gm_curr_write_mol_int_type = int_type
        temp_f = getattr(self,getattr(self, '_write_gm_' + temp_dir))
        txt2write = temp_f(**kwargs)
        if txt2write:
            self._write_directive_line(gs, temp_dir, **kwargs)
            gs.write(txt2write)

    def _write_directive(self, gs, directive, **kwargs):
        temp_f = getattr(self, self.__write_directive_v)
        temp_f(gs, directive, **kwargs)

    def _check_flag_write(self, temp_gm_object, **kwargs):
        if kwargs.get('flag_include'):
            temp_flag = self._check_flag_include(temp_gm_object, **kwargs)
            if not temp_flag:return False
        if kwargs.get('flag_if'):
            temp_flag = self._check_flag_if(temp_gm_object, **kwargs)
            if not temp_flag:return False
#        if not kwargs.get('flag_replace_define') and isinstance(temp_gm_object, Define):
        if isinstance(temp_gm_object, Define):
            temp_flag = self._check_define_inWS(temp_gm_object, **kwargs)
            if not temp_flag:
                return False
            else:
                return True
        if kwargs.get('curr_directive'):
            temp_flag = self._check_curr_directive(temp_gm_object, **kwargs)
            if not temp_flag:return False
        if kwargs.get('flag_check_directive_inWS'):
            if not (isinstance(temp_gm_object, ParsedSegments) and temp_gm_object.seg_type in ('gm_define', 'if')):
                temp_flag = self._check_directive_inWS(temp_gm_object, **kwargs)
                if not temp_flag:return False
        return True

    def _check_flag_include(self, temp_gm_object, **kwargs):
        temp_depth = len(self._gm_writing_state) - 1
        incl_depth = kwargs.get('include_depth', 0)
        temp_writing_state = self._gm_writing_state[temp_depth]
        """
        if len(temp_gm_object.gm_ord[1]) > temp_depth:
            if temp_gm_object.gm_ord[1][temp_depth] not in temp_writing_state[1]:
                temp_writing_state[1].append(temp_gm_object.gm_ord[1][temp_depth])
                self._write_include(temp_writing_state[0], temp_gm_object.gm_ord[1][temp_depth], **kwargs)
            return False
        return True
        """
        if kwargs.get('flag_include_temp_gm_object'):
            temp_include = temp_gm_object
        else:
            temp_include = None
            if hasattr(temp_gm_object, 'gm_ord') and temp_gm_object.gm_ord[1]:
                temp_include = temp_gm_object.gm_ord[1][incl_depth]
        if temp_include:
            if temp_include not in temp_writing_state[1]:
                temp_writing_state[1].append(temp_include)
                self._write_include(temp_writing_state[0][-1], temp_include, **kwargs)
            return False
        return True

    def _check_flag_if(self, temp_gm_object, **kwargs):
        if kwargs.get('flag_if_temp_gm_object'):
            temp_if = temp_gm_object
        else:
            temp_if = None
            if hasattr(temp_gm_object, 'gm_ord') and  temp_gm_object.gm_ord[2]:
                temp_if = temp_gm_object.gm_ord[2][0]
        temp_depth = len(self._gm_writing_state) - 1
        temp_writing_state = self._gm_writing_state[temp_depth]
        if temp_if:
            if temp_if not in temp_writing_state[2]:
                temp_writing_state[2].append(temp_if)
                self._write_if(temp_writing_state[0][-1], temp_if, **kwargs)
            return False
        return True

    def _check_curr_directive(self, temp_gm_object, **kwargs):
        if not hasattr(temp_gm_object, 'gm_ord'):
            return False
        return kwargs.get('curr_directive')==temp_gm_object.gm_ord[0]

    def _check_directive_inWS(self, temp_gm_object, **kwargs):
        if hasattr(temp_gm_object, 'gm_ord'):
            for ws in self._gm_writing_state:
                if temp_gm_object.gm_ord[0] in ws[4]:
                    return False
#            self.__writing_state[-1][4].append(temp_gm_object.gm_ord[0])
        return True

    def _check_define_inWS(self, temp_gm_object, **kwargs):
        if not kwargs.get('flag_replace_define'):
            for ws in self._gm_writing_state:
                if temp_gm_object in ws[3]:
                    return False
        self._gm_writing_state[-1][3].append(temp_gm_object)
        return True

    @staticmethod
    def _get_gms_from_fpaht(fpath):
        if fpath:
            if isinstance(fpath, GromacsFile):
                gs = fpath
            else:
                gs = GromacsFile(fpath, write=True)
        else:
            gs = GromacsString()
        return gs

_GromacsWriter_defs = {}
_GromacsWriter_defs['__write_gromacs_format_v'] = '__write_gromacs_format_v1'
_GromacsWriter_defs['__write_directive_v'] = '__write_directive_v1'
_GromacsWriter_defs['__write_directive_line_v'] = '__write_directive_line_v1'
_GromacsWriter_defs['__write_define_v'] = '__write_define_v1'
_GromacsWriter_defs['__write_TITLE_v'] = '__write_TITLE_v1'

GromacsWriter._add_defaults(_GromacsWriter_defs, flag_set=True)

class gmFFParser(GromacsParser):
    def __just_for_readability(self):
        self.a_type = self.get_container('a_type') # container with atom types
        self.b_a_type = self.get_container('b_a_type') # container with bond_atom types - for opls FF
        self._segments = list()
        self.get_intDB = None # from FF
        self._set_gm_interaction_type_container = None # from FF
        self._b_a_atom_wildcard = None # FF

    def parse_ff_gm(self, parse_from, parse_from_file = True, **kwargs):
        self.parse_gm(parse_from, parse_from_file, **kwargs)

    gm_FF_Defaults = gm_FF_Defaults

    def __defaults_parser_v1(self, l, *args, **kwargs):
        temp = l.split()
        if hasattr(self,'gm_ff_defaults'):
            do_warn('gm_ff_defaults already defined, overwriting...')
        self.gm_ff_defaults = gm_FF_Defaults(temp)
        self._add_gm_ord(self.gm_ff_defaults)

    def __atomtypes_parser_v1(self, l, *args, **kwargs):
        temp = l.split()
        at_id = temp[0]
        if len(temp) == 7:
            b_at_id = temp[0]
        else:
            b_at_id = temp[1]
            temp_b_at = BondAtomType(b_at_id)
            self.get_intDB().add2container(temp_b_at, create=True, db_type=set)
        temp_at = AtomType(at_id, at_id, bond_atom_id = b_at_id, vdw = [float(temp[-2]), float(temp[-1])],
                           format_type='gm')
        temp_at.element = temp[-6]
        temp_at.m = float(temp[-5])
        temp_at.p_ch = float(temp[-4])
        temp_at.p_type = temp[-3]
        self._add_gm_ord(temp_at)
        self.get_intDB().add2container(temp_at, create = True, **kwargs)
        if 'DUM' in temp_at.name:
                self.get_intDB().DUM_type = temp_at

    def __find_b_a_type(self, temp_at_id, allow_not_found = False):
        if temp_at_id == self._b_a_atom_wildcard:
            return self._b_a_atom_wildcard
        b_a_type_cont = self.get_intDB(BondAtomType, create = False, allow_not_found = 1)
        if b_a_type_cont:
            if temp_at_id in b_a_type_cont:
                return temp_at_id
        a_type_cont = self.get_intDB(AtomType, create = False, allow_not_found = 1)
        if temp_at_id in a_type_cont:
            return a_type_cont[temp_at_id]
        if allow_not_found:
            return
        raise Exception('atom type not found: '+str(temp_at_id))

        ##################################### should be split in general and specific ones
    def __parse_interaction_type(self, l, defines, int_type = None, remove_defs=False, **kwargs):
        if int_type is None:
            int_type = self._gm_curr_ff_int_type
        na = int_type.na
        temp_l = split_gen(l)
        temp_int_atoms = []
        if int_type==DihedralType:
            for i in range(2):
                temp_at_id = next(temp_l)
                temp_int_at = self.__find_b_a_type(temp_at_id)
                temp_int_atoms.append(temp_int_at)
            temp_at_id = next(temp_l)
            temp_int_at = self.__find_b_a_type(temp_at_id, allow_not_found=True)
            if temp_int_at is None:
                temp_fnc = temp_at_id
            else:
                temp_int_atoms.append(temp_int_at)
                temp_at_id = next(temp_l)
                temp_int_at = self.__find_b_a_type(temp_at_id)
                temp_int_atoms.append(temp_int_at)
                temp_fnc = next(temp_l)
        else:
            for i in range(na):
                temp_at_id = next(temp_l)
                temp_int_at = self.__find_b_a_type(temp_at_id)
                temp_int_atoms.append(temp_int_at)
            temp_fnc = next(temp_l)
        if int_type!=cmap:
            temp_int_type = int_type(fnc_type=temp_fnc, atoms=temp_int_atoms)
        else:
            temp_grid = (int(next(temp_l)), int(next(temp_l)))
            temp_int_type = int_type(fnc_type=temp_fnc, atoms=temp_int_atoms, grid_ind = temp_grid)
        self._set_gm_interaction_type_container(temp_int_type)
        if not remove_defs:
            temp_param = next(temp_l)
            if temp_param in defines:
                temp_int_type.gmdef = temp_param
                temp_int_type.set_params(defines[temp_param].v)
                extra_params = list(temp_l)
                if extra_params:
                    do_warn('not able to read all parameters, rest:', extra_params)
            else:
                temp_int_type.add_params(temp_param)
                temp_int_type.add_params(temp_l, flag_check_len = 1)
        else:
            temp_int_type.set_params(self.exchange_defs(temp_l, defines))
        self._add_gm_ord(temp_int_type)
#        if int_type == DihedralType:
#            temp_int_type._gm_check_dih_imp()
        self.get_intDB().add2container(temp_int_type, create = True, db_type=list, **kwargs)


_gmFFParser_defs = {}
_gmFFParser_defs['_defaults_parser'] = '__defaults_parser_v1'
_gmFFParser_defs['_atomtypes_parser'] = '__atomtypes_parser_v1'

for temp_dir in direction_interaction_type_map:
    temp_dir_fnc_attr = '_'+temp_dir+'_parser'
    if temp_dir_fnc_attr not in _gmFFParser_defs:
        _gmFFParser_defs[temp_dir_fnc_attr] = '__parse_interaction_type'

#_gmFFParser_defs['_cmaptypes_parser'] = '__cmaptypes_parser_v1'
#_gmFFParser_defs['_nonbond_params_parser'] = '__nonbond_params_parser_v1'
#_gmFFParser_defs['_pairtypes_parser'] = '__pairtypes_parser_v1'

gmFFParser._add_defaults(_gmFFParser_defs, flag_set=1)
fnc_default = {'gb': '2', 'ga': '2', 'gd': '1', 'gi': '2'}


class gmFFWriter(GromacsWriter):

    __a_type_format1 = ' {:>4} {:>4} {:>10} {:>10} {:>5} {:12.7e} {:12.7e}\n'
    __a_type_format2 = ' {:>9} {:>3} {:>3} {:>10} {:>10} {:>5} {:12.7e} {:12.7e}\n'

    __int_type_atom_format = '{:<6} '
    __int_type_fnc_format = '{:<1}'
    __int_type_p_format = {str:' {:12s}', float:' {:12.7e}', int: ' {:12d}'}
    __int_type_gmdef_format = ' {:>6}'

    def _write_gm_ff_defines(self, gs, **kwargs):
        defs2write = []
        for define in getattr(self, 'defines', []):
            if define.startswith('_FF_'):
                defs2write.append(define)
        kwargs['flag_replace_define'] = False
        self._write_define(gs, defs = defs2write, **kwargs)

    def __write_gm_defaults_v1(self, **kwargs):
        if kwargs.get('flag_ff_defines'):
            gs = self._gm_writing_state[-1][0][-1]
            self._write_gm_ff_defines(gs, **kwargs)
        temp_txt = ''
        if self._check_flag_write(self.gm_ff_defaults, **kwargs):
            temp_txt = self.gm_ff_defaults.write_defs()
        return temp_txt

    def __write_gm_atomtypes_v1(self, **kwargs):
        """
        :param kwargs:
        :return:
        """
        txt2write = ''
        a_type_cont = self.get_intDB(AtomType)
        b_a_type_cont = self.get_intDB(BondAtomType, create=False, allow_not_found=1)
        for at_id in a_type_cont:
            temp_a_type = a_type_cont[at_id]
            flag_write = self._check_flag_write(temp_a_type, **kwargs)
            if not flag_write:continue
            temp_v = []
            temp_v.append(get_id_string(temp_a_type, flag_id='gm_id'))
            if b_a_type_cont:
                temp_v.append(temp_a_type.b_id)
            temp_v.append(temp_a_type.element)
            temp_v.append(temp_a_type.m)
            temp_v.append(temp_a_type.p_ch)
            temp_v.append(temp_a_type.p_type)
            temp_v.extend(temp_a_type.vdw)
            if b_a_type_cont:
                txt = self.__a_type_format2.format(*temp_v)
            else:
                txt = self.__a_type_format1.format(*temp_v)
            txt2write += txt
        return txt2write

    def _check_vdw_pair_gm_params(self):
        ff = self.get_intDB()
        if not hasattr(ff,'vdw'):
            self.generate_vdw()
        for ii, at_pair in enumerate(ff.vdw):
            temp_vdwpair_vdw = ff.vdw[at_pair]
            temp_vdwtype2 = self.find_interaction_type(at_pair, vdWType, self.gm_ff_defaults[0], use_vdw = False,
                                                       allow_not_found = 1, create = False)
            temp_pairtype2 = self.find_interaction_type(at_pair, PairType, self.gm_ff_defaults[0], use_vdw = False,
                                                       allow_not_found = 1, create = False)
            flag_generated = [False, False]
            if temp_vdwtype2:
                temp_vdwtype2 = temp_vdwtype2[0]
            if temp_pairtype2:
                temp_pairtype2 = temp_pairtype2[0]
            if not temp_vdwtype2 or (not temp_pairtype2 and temp_vdwpair_vdw[1]):
                temp_vdwpair_type2 = self._generate_vdw_atpair_gm(at_pair[0].id, at_pair[1].id, replace=-1)
                if not temp_vdwtype2:
                    temp_vdwtype2 = temp_vdwpair_type2[0]
                    flag_generated[0] = True
                if not temp_pairtype2 and temp_vdwpair_vdw[1]:
                    temp_pairtype2 = temp_vdwpair_type2[1]
                    flag_generated[1] = True
            for i in range(2):
                temp_new_vdw_pair_type = (temp_vdwtype2, temp_pairtype2)[i]
                if temp_vdwpair_vdw[i]:
                    if not temp_new_vdw_pair_type:
                        self._set_gm_interaction_type_container(temp_vdwpair_vdw[i])
                        ff.add2container(temp_vdwpair_vdw[i], replace=1, create=True, db_type=list)
                    elif not temp_vdwpair_vdw[i].check_eq_params(temp_new_vdw_pair_type.p):
                        if flag_generated[i]:
                            self._set_gm_interaction_type_container(temp_vdwpair_vdw[i])
                            ff.add2container(temp_vdwpair_vdw[i], replace=1, create=True, db_type=list)
                        else:
                            temp_new_vdw_pair_type.p = temp_vdwpair_vdw[i].p

    def __write_interaction_type_v1(self, **kwargs):
        int_type = self._gm_curr_write_ff_int_type
        temp_cont = self.get_intDB(int_type, cont_pref = self._gm_cont_pref, create = False, allow_not_found = 1)
        if not temp_cont:return
        txt2write = ''
        for temp_int_type in temp_cont:
            flag_write = self._check_flag_write(temp_int_type, **kwargs)
            if not flag_write:continue
            defines = None
            if kwargs.get('flag_use_define'):
                defines = self.defines
            temp_txt = temp_int_type._write_gm_int_type(defines = defines, **kwargs)
            txt2write += temp_txt + '\n'
        return txt2write

    def write_ff_itp(self, f_path = None, **kwargs):
        """
        :param f_path:
        :param kwargs:
             flag_split_non_bonded = False
             sep2include_files - customizable
             flag_close = True (close gromacs stream)
        :return:
        """
        gs = self._get_gms_from_fpaht(f_path)
        if kwargs.get('flag_split_non_bonded'):
            kwargs['sep2include_files'] = [(None, 1), ('ffnonbonded.itp',2), ('ffbonded.itp',1)]
        self.write_gromacs_format(gs, ff_directives, **kwargs)
        if f_path and kwargs.get('flag_close', True):
            gs.f.close()


_gmFFWriter_defs = {}
_gmFFWriter_defs['_write_gm_defaults'] = '__write_gm_defaults_v1'
_gmFFWriter_defs['_write_gm_atomtypes'] = '__write_gm_atomtypes_v1'
for temp_dir in direction_interaction_type_map:
    _gmFFWriter_defs['_write_gm_' + temp_dir] = '__write_interaction_type_v1'
#_gmFFWriter_defs['_write_gm_cmaptypes'] = '__write_cmaptypes_v1'

gmFFWriter._add_defaults(_gmFFWriter_defs, flag_set=1)


"""
read / write for rtp???

"""

class gmFragmentMoleculeIO(GromacsParser, GromacsWriter):

    def __just_for_readability(self):
        self.ff = None # FF
        #self.find_interaction_type = None # from FF
        #self.get_intDB = None # from FF
        self.add_residue = None # Molecule
        self.Atom = None # Molecule
        self.ChargeGroup = None # Molecule
        self.cg = None # Molecule
        self.atoms = None # Molecule

    @staticmethod
    def __find_a_type_gm(ff, atom_type_str):
        for a_type in ff.a_type.values():
            a_type_gm_id = getattr(a_type, 'gm_id', None)
            if a_type_gm_id and a_type_gm_id==atom_type_str:
                return a_type

    def __atoms_parser_v1(self, l, defines, **kwargs):
        temp = split_gen(l)
        atom_number = int(next(temp))
        atom_type_str = next(temp)
        ff = self.ff.get_intDB()
        atom_type = ff.get_item(atom_type_str, AtomType, allow_not_found = 1)
        if atom_type is None:
            atom_type = self.__find_a_type_gm(ff, atom_type_str)
        if atom_type is None:
            atom_type = atom_type_str
        res_n = int(next(temp))
        res_name = next(temp)
        temp_res = self.add_residue(res_n, res_name)
        atom_name = next(temp)
        cg_n = int(next(temp))
        if hasattr(self, 'cg') and self.cg:
            if self.cg[-1].n != cg_n:
                temp_cg = self.get_item(None, self.ChargeGroup, create=True, create_container=True, db_type = list)
                temp_cg.n = cg_n
        else:
            temp_cg = self.get_item(None, self.ChargeGroup, create=True, create_container=True, db_type=list)
            temp_cg.n = cg_n
        p_ch = float(next(temp))
        try:
            mass = float(next(temp))
        except:
            mass = atom_type.m
        temp_atom = self.Atom(atom_number, name=atom_name)
        temp_atom.a_type = atom_type
        temp_atom.p_ch = p_ch
        temp_atom.m = mass
        temp_atom._generate_self_state()
        temp_res.gm_id = temp_res.id
        temp_res.add_atom(temp_atom)
        try:
            atom_type_B_str = next(temp)
            atom_type_B = self.ff.get_intDB().get_item(atom_type_B_str, AtomType, allow_not_found = 1)
            if atom_type_B is None:
                atom_type_B = atom_type_B_str
            partial_charge_B = float(next(temp))
            mass_B = float(next(temp))
            temp_atom.add_atom_state(atom_type_B, partial_charge_B, mass_B)
        except:
            pass
        self.cg[-1].add_atom(temp_atom)
        self._add_gm_ord(temp_atom, **kwargs)
        temp_atom.gm_id = temp_atom.id
        self.add2container(temp_atom, create=True)

    def __parse_atoms_gm_interaction(self, l, interaction_type):
        na = interaction_type.na
        temp_l = split_gen(l)
        temp_int_atoms = []
        for i in range(na):
            temp_at_id = int(next(temp_l))
            temp_int_at = self.atoms[temp_at_id]
            temp_int_atoms.append(temp_int_at)
        temp = list(temp_l)
        return temp_int_atoms, temp, na

    def _parse_gm_interaction(self, l, defines, interaction_type = None, default_fnc_type=False, **kwargs): # not the nicest solution
        if interaction_type is None:
            interaction_type = self._gm_curr_int_type
        temp_int_atoms, temp, na = self.__parse_atoms_gm_interaction(l, interaction_type)
        if not default_fnc_type:
            if len(temp) < 1:
                raise Exception("missing function type while reading ", interaction_type.__name__)
            fnc_type = temp.pop(0)
        else:
            fnc_type = default_fnc_type
        temp_interaction_type = interaction_type(fnc_type = fnc_type).__class__
        temp_interaction = Interaction(temp_interaction_type, atoms=temp_int_atoms)
        if not temp:
            temp_int_atom_types = [temp_at.a_type for temp_at in temp_int_atoms]
            temp_state = self.ff.find_interaction_type(temp_int_atom_types, temp_interaction_type, fnc_type, **kwargs)
            temp_interaction.add_state(temp_state[0])
        else:
            if temp_interaction_type!=VirtualSitenType:
                while temp:
                    temp_interaction.add_state(fnc_type = fnc_type)
                    temp_n_param = len(temp_interaction.states[-1].p_type)
                    if temp[0] in defines:
                        temp_def = temp.pop(0)
                        temp_params = defines[temp_def].v
                        temp_interaction.states[-1].add_params(temp_params)
                        temp_interaction.states[-1].gmdef = temp_def
                    else:
                        temp_params, temp = temp[:temp_n_param], temp[temp_n_param:]
                        temp_interaction.states[-1].add_params(temp_params)
            else:
                temp_interaction.w = []
                if fnc_type != '3':
                    for temp_at_id in temp:
                        additional_atom = self.atoms[temp_at_id]
                        temp_interaction.atoms.append(additional_atom)
                        if fnc_type=='1':
                            temp_interaction.w.append(1)
                        elif fnc_type=='2':
                            temp_interaction.w.append(additional_atom.m)
                        else:
                            raise Exception('fnc_type not known', fnc_type)
                else:
                    if len(temp)%2:
                        raise Exception('number of paramerets wrong: len(temp) % 2 should be 0 ', len(temp))
                    for i in range(int(len(temp)/2)):
                        at, temp_w = self.atoms[temp[i * 2]], float(temp[i * 2 + 1])
                        temp_interaction.atoms.append(at)
                        temp_interaction.w.append(temp_w)
        self._add_gm_ord(temp_interaction, **kwargs)
        self.add2container(temp_interaction, create = True, db_type=list)

    def __exclusions_parser_v1(self, l, defines, **kwargs):
        interaction_type = kwargs.get('interaction_type', ExclusionType)
        temp_l = split_gen(l)
        temp_at_id = int(next(temp_l))
        at1 = self.atoms[temp_at_id]
        for temp_at_id in temp_l:
            at2 = self.atoms[int(temp_at_id)]
            #at1.add_excl(self.atoms[int(temp_at_id)])
            temp_interaction = Interaction(interaction_type, atoms=(at1, at2))
            temp_interaction.add_state(fnc_type = None, params = (True,))
            self._add_gm_ord(temp_interaction, **kwargs)
            self.add2container(temp_interaction, create=True, db_type=list)
            # exclusion pair type
            temp_interaction = Interaction(ExclusionPairType, atoms=(at1, at2))
            temp_interaction.add_state(fnc_type = None, params = (True,))
            self._add_gm_ord(temp_interaction, **kwargs)
            self.add2container(temp_interaction, create=True, item_id=frozenset((at1, at2)), replace = -1)
            self.add_atom_pair2EP_l(at1, at2)

    def __virtualsitesn_parser_v1(self, l, defines, **kwargs):
         # it's in _parse_gm_interaction, but should be separated...
        pass

#writing part
######################################################################## PTP missing!!!!!

    def write_itp(self, f_path, directives = None, **kwargs):
        gs = self._get_gms_from_fpaht(f_path)
        self.write_molecule_type(gs, **kwargs)
        if f_path and kwargs.get('flag_close', True):
            gs.f.close()

    def write_molecule_type(self, gs, directives = None, **kwargs):
        """
        :param gs:
        :param directives:
        :param kwargs:
             flag_segment_order = True
             flag_general_defines = False(from top e.g) - generates defines2write
             defines2write - container of defines
        :return:
        """
        if self.id is None:
            self.id = kwargs.get('mol_id', self.id_pref + '1')
        if kwargs.get('flag_generate_excl_pairs'):
            self._gm_generate_excl_pairs(**kwargs)
        kwargs = dict(kwargs)
        top_state, other_state = self._get_top_other_state(**kwargs)
        top_other_state = dict(top_state=top_state, other_state=other_state)
        kwargs.update(top_other_state)
        kwargs['flag_ff_defines'] = 0
        if kwargs.get('flag_segment_order', True):
            kwargs['gm_segments_writing_stack'] = list()
            if hasattr(self.ff, '_segments') and hasattr(self, 'gm_ord'):
                if self.gm_ord[0] in self.ff._segments:
                    pos = self.ff._segments.index(self.gm_ord[0])
                    kwargs['gm_segments_writing_stack'] = list(self.ff._segments[pos:])
            if hasattr(self, '_segments'):
                kwargs['gm_segments_writing_stack'].extend(list(self._segments))
        if kwargs.get('flag_general_defines', False):
            if not kwargs.get('defines2write'):
                kwargs['defines2write'] = getattr(self.ff, 'defines')
        if not kwargs.get('gm_segments_writing_stack'):
            kwargs['flag_defines'] = 0
#        kwargs['check_self_flag_write'] = 1
#        if self._check_flag_write(self, **kwargs):
        if directives is None:
            directives = mol_directives
        self.write_gromacs_format(gs, directives, **kwargs)

    def __write_moleculetype_v1(self, **kwargs):
        if self._check_flag_write(self, **kwargs):
            return '{:<20} {:}\n'.format(self.id, self.nrexcl)
        return ''

    __atom_format = '{:>5}'
    __atom_type_format = '{:>6s}'
    __atoms_line_format1 = '{:>6}{:>11}{:>7}{:>7}{:>7}{:>7}'
    __pch_format = ' {:>10f}'
    __m_format = ' {:>10f}'
    __fnc_type_format = '{:>5s}'

    def __write_atoms_v1(self, **kwargs):
        txt2write = ''
        a_type_kwargs, m_kwargs, pch_kwargs = self._get_atom_kwargs(**kwargs)
        ptp_atom_text = ' ' + self.__atom_type_format + ' ' + self.__pch_format + ' ' + self.__m_format
        for at_id in self.atoms:
            at = self.atoms[at_id]
            flag_write = self._check_flag_write(at, **kwargs)
            if not flag_write:continue
            a_type_ptp = self.get_ptp_states(at, **a_type_kwargs)
            m_ptp = self.get_ptp_states(at, **m_kwargs)
            p_ch_ptp = self.get_ptp_states(at, **pch_kwargs)
            a_type_ptp, m_ptp, p_ch_ptp, flag_ptp = self.get_atom_ptp_states(at, a_type_ptp, m_ptp, p_ch_ptp, **kwargs)
            txt2write += self.__atoms_line_format1.format(get_id_string(at, flag_id='gm_id'),
                                get_id_string(a_type_ptp[0], flag_id='gm_id'), get_id_string(at.res, flag_id='gm_id'),
                                                          at.res.name, at.name, at.cg.n)
            txt2write += (' ' + self.__pch_format + ' ' + self.__m_format).format(p_ch_ptp[0], m_ptp[0])
            ## PTP
            if flag_ptp:
                txt2write += ptp_atom_text.format(get_id_string(a_type_ptp[1], flag_id='gm_id'), p_ch_ptp[1], m_ptp[1])
            txt2write += '\n'
        return txt2write

    def __write_interactions_v1(self, **kwargs):
        int_type = kwargs.get('interaction_type', self._gm_curr_write_mol_int_type)
        int_container = self.get_container(int_type.find_int_container2write(), create = False, allow_not_found = 1)
        if not int_container:return
        txt2write = ''
        for temp_interaction in int_container:
            flag_write = self._check_flag_write(temp_interaction, **kwargs)
            if not flag_write:continue
            temp_txt = (self.__atom_format + ' ') * len(temp_interaction.atoms)
            txt2write += temp_txt.format(*[get_id_string(at, flag_id='gm_id') for at in temp_interaction.atoms])
            temp_int_ptp = self.get_ptp_states(temp_interaction)
            if temp_int_ptp:
                state4writing, ptp_state4writing = temp_int_ptp
            else:
                ptp_state4writing = None
                state4writing = self._get_state(temp_interaction, **kwargs)
            """
            if hasattr(temp_interaction, 'state'):
                state4writing = temp_interaction.state
            else:
                state4writing = self.get_state(temp_interaction.states, **kwargs)
            """
            if hasattr(state4writing, 'gm_int_type'):
                temp_gm_int_type = state4writing.gm_int_type
            else:
                temp_gm_int_type = state4writing
            if ptp_state4writing and hasattr(ptp_state4writing, 'gm_int_type'):
                temp_ptp_gm_int_type = ptp_state4writing.gm_int_type
            else:
                temp_ptp_gm_int_type = ptp_state4writing
            """
            if temp_interaction.states[0].fnc_type=='gr_fnc':
                temp_gm_int_type = getattr(temp_interaction.states[0],'gm_int_type', None)
                if temp_gm_int_type is None:
                    temp_gm_int_type = temp_interaction.states[0].gr2gm(**kwargs)
            else:
                temp_gm_int_type = temp_interaction.states[0]
            """
            if temp_gm_int_type.fnc_type:
                txt2write += self.__fnc_type_format.format(temp_gm_int_type.fnc_type)
            if kwargs.get('write_params', True):
                flag_write_params = 1
                if not kwargs.get('write_params_explicitly'):
                    temp_int_type = self.ff.find_interaction_type(temp_interaction.atoms, temp_interaction.int_type,
                        temp_gm_int_type.fnc_type, create = False, allow_not_found=1, allow_multiple_matches = True)
                    """
                    temp_int_type = self.ff.find_interaction_type(temp_interaction.atoms, temp_interaction.int_type,
                        temp_gm_int_type.fnc_type, int_type_container = [temp_interaction._sol[0]],
                                                                  create = False, allow_not_found=1)
                    """
                    if temp_int_type and temp_int_type[0].check_eq_params(temp_gm_int_type.p):
                        flag_write_params = 0
                if flag_write_params:
                    defines = None
                    if kwargs.get('flag_use_define'):
                        defines = kwargs.get('defines', getattr(self.ff, 'defines', None))
                    txt2write += ' ' + temp_gm_int_type._write_gm_params(defines = defines, **kwargs)
                # PTP interaction
                if temp_ptp_gm_int_type:
                    flag_write_params = 1
                    if not kwargs.get('write_params_explicitly'): ######################### this part has to be adjusted for PTP atom types
                        temp_int_type = self.ff.find_interaction_type(temp_interaction.atoms, temp_interaction.int_type,
                        temp_ptp_gm_int_type.fnc_type, create = False, allow_not_found=1, allow_multiple_matches = True)
                        if temp_int_type and temp_int_type[0].check_eq_params(temp_ptp_gm_int_type.p):
                            flag_write_params = 0
                    if flag_write_params:
                        defines = None
                        if kwargs.get('flag_use_define'):
                            defines = kwargs.get('defines', getattr(self.ff, 'defines', None))
                        txt2write += ' ' + temp_ptp_gm_int_type._write_gm_params(defines = defines, **kwargs)
            txt2write += '\n'
        return txt2write

    def __write_exclusions_v1(self, **kwargs):
        kwargs['write_params'] = False
        kwargs['interaction_type']= ExclusionType
        temp = self.__write_interactions_v1(**kwargs)
        return self.__write_interactions_v1(**kwargs)


_gmFragmentMoleculeIO_defs = {}
_gmFragmentMoleculeIO_defs['_atoms_parser'] = '__atoms_parser_v1'
_gmFragmentMoleculeIO_defs['_exclusions_parser'] = '__exclusions_parser_v1'
for temp_dir in direction_interaction_map:
    temp_dir_fnc_attr = '_'+temp_dir+'_parser'
    if temp_dir_fnc_attr not in _gmFragmentMoleculeIO_defs:
        _gmFragmentMoleculeIO_defs[temp_dir_fnc_attr] = '_parse_gm_interaction'

mol_directives = [['moleculetype'], ['atoms', 'bonds', 'exclusions']]
for temp_dir in direction_interaction_map:
    if temp_dir not in mol_directives[-1]:
        mol_directives[-1].append(temp_dir)

_gmFragmentMoleculeIO_defs['_write_gm_moleculetype'] = '__write_moleculetype_v1'
_gmFragmentMoleculeIO_defs['_write_gm_atoms'] = '__write_atoms_v1'
_gmFragmentMoleculeIO_defs['_write_gm_exclusions'] = '__write_exclusions_v1'

for temp_dir in direction_interaction_map:
    _gmFragmentMoleculeIO_defs['_write_gm_' + temp_dir] = '__write_interactions_v1'

gmFragmentMoleculeIO._add_defaults(_gmFragmentMoleculeIO_defs, flag_set=1)


class gmTopologyIO(GromacsParser, GromacsWriter):
    def __just_for_readability(self):
        self.Molecule = None
        self.MoleculeType = None
        self.molecule_types = None

    def parse_top_gm(self, parse_from, parse_from_file = True, **kwargs):
        self.parse_gm(parse_from, parse_from_file, **kwargs)
        for mt in self.molecule_types.values():
            mt.add_exclusions_neigh(nexcl = mt.nrexcl)
            pairs_cont = mt.get_container(PairType, flag_class = True, allow_not_found = True)
            if pairs_cont:
                for temp_pair in pairs_cont:
                    temp_interaction = Interaction(ExclusionPairType, atoms=temp_pair.atoms)
                    temp_interaction.add_state(fnc_type = None, params = (2,))
                    mt.add2container(temp_interaction, create=True, item_id=frozenset(temp_pair.atoms), replace = True)
                    at1, at2 = temp_pair.atoms
                    mt.add_atom_pair2EP_l(at1, at2)

    def __moleculetype_parser_v1(self, l, defines, **kwargs):
        temp_l = split_gen(l)
        temp_mol = self.MoleculeType(next(temp_l), int(next(temp_l)), ff = self)
        self._add_gm_ord(temp_mol)
        self.add2container(temp_mol, create=True)
        self._curr_bbmol = temp_mol

    def __system_parser_v1(self, l, defines, **kwargs):
        curr_directive = self.gm_curr_stack[0]
        file_type = kwargs.get('file_type')
        if hasattr(curr_directive, DescriptionPart.container2write):
            temp_sys_title = getattr(curr_directive, DescriptionPart.container2write)
        else:
            temp_sys_title = DescriptionPart(form='gm', file_type=file_type)
            setattr(curr_directive, DescriptionPart.container2write, temp_sys_title)
        temp_sys_title.lines.append(l)
        temp_cont = self.get_container(temp_sys_title, create=True, db_type=list, flag_item=True)
        self._add_gm_ord(temp_sys_title)
        temp_sys_title.flag_gm_system = True
        self.add2container(temp_sys_title, create=True, db_type=list, replace = -1)

    def __molecules_parser_v1(self, l, defines, **kwargs):
        temp_l = split_gen(l)
        temp_mol_type = self.get_item(next(temp_l), self.MoleculeType)
        temp_mol = self.Molecule(temp_mol_type, int(next(temp_l)))
        self._add_gm_ord(temp_mol)
        self.add2container(temp_mol, create=True, db_type=list)

#writing
    __molecule_name_format = '{:<20}'

    #def __write_gm_moleculetype_v1(self, gs, directives = None, **kwargs):
    def __write_gm_moleculetype_v1(self, **kwargs):
        kwargs['gm_writing_state'] = self._gm_writing_state
        molecules2write = kwargs.get('molecules2write', None)
        if molecules2write is None:
            molecules2write = self.molecule_types.keys()
        sep_mol2itp = kwargs.get('sep_mol2itp', {})
        assert isinstance(sep_mol2itp, dict) or sep_mol2itp in (0, 1, False, True)
        if sep_mol2itp is not None and not isinstance(sep_mol2itp, dict) and sep_mol2itp:
            sep_mol2itp = {}
            sep_mol2exclude = kwargs.get('sep_mol2exclude', set())
            if kwargs.get('flag_sep_mol_exclude_included', True):
                for mol_type_id in molecules2write:
                    mol_type = self.molecule_types[mol_type_id]
                    if hasattr(mol_type, 'gm_ord') and mol_type.gm_ord[1]:
                        sep_mol2exclude.add(mol_type_id)
            for mol_type_id in molecules2write:
                if mol_type_id not in sep_mol2exclude:
                    sep_mol2itp[mol_type_id] = sep_mol2itp.get(mol_type_id, mol_type_id + '.itp')
        for mol_type_id in molecules2write:
            mol_type = self.molecule_types[mol_type_id]
            if mol_type_id in sep_mol2itp:
                temp_fpath = sep_mol2itp[mol_type_id]
                new_gs = self._add_include_gs2writing_state(temp_fpath, **kwargs)
#                self._write_include(self._gm_writing_state[-1][0][-1], temp_fpath, **kwargs)
#                new_gs = GromacsFile(temp_fpath, True)
#                self._add_gs2writing_state(new_gs)
            mol_type.write_molecule_type(self._gm_writing_state[-1][0][-1], **kwargs)
            if mol_type_id in sep_mol2itp:
                gs = self._remove_gs_from_writing_state()
#                self._gm_writing_state[-1][0].pop(-1)

    def _write_gm_top_moleculetype(self, **kwargs):
        temp_f = getattr(self, self.__write_gm_moleculetype_v)
        temp_f(**kwargs)

    def __get_gm_system_txt(self):
        return 'molecular system'

    def __write_system_v1(self, **kwargs):
        txt2write = ''
        temp_cont = self.get_container(DescriptionPart, flag_class=True, allow_not_found = 1)
        for dp in temp_cont:
            flag_dp = 0
            if getattr(dp, 'flag_gm_system', None):
                flag_dp = 1
            else:
                temp_gm_ord = getattr(dp, 'gm_ord', None)
                if temp_gm_ord  and temp_gm_ord[0].seg_type == 'gm_directive' and temp_gm_ord[0].seg_name == 'system':
                    flag_dp = 1
            if flag_dp and dp.lines:
                txt2write+=dp.lines[0]
        if not txt2write:
            txt2write = self.__get_gm_system_txt()
        return txt2write

    def __write_molecules_v1(self, **kwargs):
        txt2write = ''
        for m in self.molecules:
            if m.n:
                txt2write += (self.__molecule_name_format + ' {:}\n').format(m.mol_type.id, m.n)
        return txt2write

    def write_top_gm(self, f_path = None, **kwargs):
        """
        :param f_path:
        :param kwargs:
            see write_ff_itp
            sep_ff2itp - separate ff into itp file(s)
            sep_mol2itp - separate molecules to itp files
                sep_mol2itp = {mol_id1:f_path1}
        :return:
        """
        gs = self._get_gms_from_fpaht(f_path)
        ff_kwargs = copy.deepcopy(kwargs)
        ff_kwargs['flag_del_writing_state'] = False
        ff_kwargs['flag_close'] = False
        sep_ff2itp = kwargs.get('sep_ff2itp')
        if sep_ff2itp:
            if sep_ff2itp == 1:
                sep_ff2itp = 'ff.itp'
            assert isinstance(sep_ff2itp, str)
            ff_kwargs['last_inc_file'] = sep_ff2itp
        self.write_ff_itp(gs, **ff_kwargs)
#        self.write_gromacs_format(gs, ff_directives, **kwargs)
        self._write_gm_top_moleculetype(**kwargs)
        kwargs['gm_writing_state'] = self._gm_writing_state
        kwargs['gm_segments_writing_stack'] = self._segments_writing_stack
        kwargs['flag_title'] = 0
        self.write_gromacs_format(gs, top_directives, **kwargs)
        if f_path and kwargs.get('flag_close', True):
            gs.f.close()

_gmTopologyIO_defs = {}
_gmTopologyIO_defs['_moleculetype_parser'] = '__moleculetype_parser_v1'
_gmTopologyIO_defs['_molecules_parser'] = '__molecules_parser_v1'
_gmTopologyIO_defs['_system_parser'] = '__system_parser_v1'

_gmTopologyIO_defs['_write_gm_system'] = '__write_system_v1'
_gmTopologyIO_defs['__write_gm_moleculetype_v'] = '__write_gm_moleculetype_v1'
_gmTopologyIO_defs['_write_gm_molecules'] = '__write_molecules_v1'

gmTopologyIO._add_defaults(_gmTopologyIO_defs, flag_set=1)






















