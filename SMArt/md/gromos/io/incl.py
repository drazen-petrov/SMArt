from itertools import chain
from SMArt.incl import os, copy, np, do_warn
from SMArt.incl import Defaults, OrderedDict, GeneralContainer, FileStream, StringStream, split_gen, get_id_string
from SMArt.md.incl import *
from SMArt.md.incl import _gr_title_gm_system, _file_type_ext_map, _add_DescriptionPart
from ..incl import GromosDefaults, GromosBlockNames

#### GENERAL CLASSES / FUNCTIONS

class GromosStream(Defaults):
    _comm_indicators = ('#',)
    _blockend = 'END'
    __buffer_l = ''

    def __just_for_readability(self):
        self.f = None
        self.lines = None #this is implemented in GromosFile or GromosString

    def block_lines(self):
        """generator yielding line by line of a block (comments removed)"""
        for i in self.lines():
            if not i.startswith(self._blockend):
                yield i
            else:
                self.__buffer_l = ''
                return

    @property
    def block_split_fnc(self):
        """generator yielding elements of a block split by empty spaces (' ' or '\t' ...)"""
        if self.__buffer_l == '':
            self.__buffer_l = split_gen('')
        while 1:
            try:
                yield next(self.__buffer_l)
            except StopIteration:
                try:
                    self.__buffer_l = split_gen(next(self.block_lines()))
                except StopIteration:
                    self.__buffer_l = ''
                    break

    @property
    def report_nl_block_split_fnc(self):
        if self.__buffer_l == '':
            self.__buffer_l = split_gen('')
        try:
            temp = next(self.__buffer_l)
            self.__buffer_l = chain(split_gen(temp), self.__buffer_l)
            return False
        except StopIteration:
            return True

    def end_block_split(self, bl_name=''):
        """reads the last END of the block - basically a workaround"""
        temp = list(self.block_split_fnc)
        self.__buffer_l = ''
        if temp:
            raise Exception("more arguments than expected in block: " + bl_name)

    def next_block(self):
        """goes to the next block - reads lines till the end of the block"""
        for i in self.lines():
            if i.startswith(self._blockend):
                return

    @property
    def block_names(self):
        """reads the whole file and yields block names"""
        yield (next(self.lines()))
        for i in self.lines():
            if i.startswith(self._blockend):
                try:
                    yield (next(self.lines()))
                except StopIteration:
                    pass

# writing part
    def write_block(self, bl_name, text2write):
        """text2write: str or list of lines"""
        s2w = ''
        s2w+=bl_name + '\n'
        #self.f.write(bl_name + '\n')
        try:## python2 unicode
            if type(text2write) == unicode:
                text2write = str(text2write)
        except:
            pass
        if type(text2write) == str:
            if not text2write.endswith("\n"):
                text2write += "\n"
            s2w += text2write
            #self.f.write(text2write)
        elif type(text2write) == list:
            for i in text2write:
                if not i.endswith("\n"):
                    i += "\n"
                s2w += i
                #self.f.write(i)
        else:
            raise TypeError("text2write expected to be str or list - provided " + type(text2write).__name__)
        s2w += self._blockend + '\n'
        #self.f.write(self._blockend + '\n')
        self.f.write(s2w)
        if isinstance(self, GromosString):
            self.s+=s2w

GromosStream._add_defaults({'_blockend':'END', '_comm_indicators':('#',)})


class GromosFile(FileStream, GromosStream):
    """handler class for gromos files - reading / writing functions"""
    pass


class GromosString(StringStream, GromosStream):
    """handler class for gromos string - reading / writing functions"""
    pass


# parsing & writing classes and functions
# these classes assume that they would be used as subclasses, so some functions are not explicitly available here...

##### General parsing

class GromosParser(GeneralContainer, GromosDefaults):
    """contains basic functions that allow reading of different blocks"""
    @staticmethod
    def __parse_from_source(parsing_from, parse_from_file=True, **kwargs):
        if isinstance(parsing_from, (GromosFile, GromosString)):
            return parsing_from
        if parse_from_file:
            return GromosFile(parsing_from, **kwargs)
        return GromosString(parsing_from)

    def parse_gro(self, parse_from, parse_from_file=True, **kwargs):
        """
        parsing gromos format
        :param parse_from:
        :param parse_from_file:
        :param kwargs:

        :return:
        """
        temp_f = getattr(self, self.__parse_gro_v)
        temp_f(parse_from, parse_from_file, **kwargs)

    def __find_parse_fnc(self, bl_name):
        if hasattr(self, '_' + bl_name + '_parser'):
            fnc_name = getattr(self, '_' + bl_name + '_parser')
            return getattr(self, fnc_name)
        else:
            do_warn('unknown block read in as text: ' + bl_name)
            return self.__read_unknown_block

    def __parse_gro_v1(self, parsing_from, parse_from_file, **kwargs):
        gs = self.__parse_from_source(parsing_from, parse_from_file, **kwargs)
        if not hasattr(self, 'parsed_from'):
            self.parsed_from = []
        self.parsed_from.append(gs)
        if not hasattr(self, 'undefined_bl'):
            self.undefined_bl = OrderedDict()
        for i in gs.lines():
            bl_name = i.strip()
            temp_fnc = self.__find_parse_fnc(bl_name)
            temp_fnc(gs, bl_name=bl_name, **kwargs)
        gs._remove_f()

    def __read_unknown_block(self, gs, bl_name, **kwargs):
        self.undefined_bl[bl_name] = list(gs.block_lines())

    def __TITLE_v1(self, gs, *args, **kwargs):
        temp_title = DescriptionPart(form = 'gr', **kwargs)
        if kwargs.get('flag_add_source', True):
            temp_title._add_source(gs.f_path, additional_txt = 'GROMOS')
        temp_title.lines.extend(list(gs.block_lines()))
        self.add2container(temp_title, create=True, db_type=list)

    def add_gr_title(self, title_lines, **kwargs):
        temp_title = DescriptionPart(form = 'gr', **kwargs)
        temp_title.lines.extend(title_lines)
        self.add2container(temp_title, create=True, db_type=list)

    def __PHYSICALCONSTANTS_v1(self, gs, *args, **kwargs):
        setattr(self, getattr(self._gr_block_names, kwargs['bl_name']), list(gs.block_lines()))

    def __oneliner_block_v1(self, gs, *args, **kwargs):
        temp = list(gs.block_lines())
        if len(temp) != 1:
            raise Exception('expected 1 line block and got ' + str(len(temp)) + 'lines')
        setattr(self, getattr(self._gr_block_names, kwargs['bl_name']), temp[0].strip())

_GromosParser_defs = {}
_GromosParser_defs['__parse_gro_v'] = '__parse_gro_v1' # this will take __parse_gro_v1 function
_GromosParser_defs['_TITLE_parser'] = '__TITLE_v1'
_GromosParser_defs['_FORCEFIELD_parser'] = '__oneliner_block_v1'
_GromosParser_defs['_MAKETOPVERSION_parser'] = '__oneliner_block_v1'
_GromosParser_defs['_PHYSICALCONSTANTS_parser'] = '__PHYSICALCONSTANTS_v1'
_GromosParser_defs['_oneliner_block'] = '__oneliner_block_v1'

GromosParser._add_defaults(_GromosParser_defs, flag_set=1)

##### General writing

class GromosWriter(GeneralContainer, GromosDefaults):
    """contains basic functions that allow writing of different blocks"""

    @staticmethod
    def _get_grs_from_fpaht(fpath):
        if fpath:
            if isinstance(fpath, GromosFile):
                gs = fpath
            else:
                gs = GromosFile(fpath, write=True)
        else:
            gs = GromosString()
        return gs

    def __write_unknown_blocks(self, gs, *blocks, **kwargs):
        for bl in blocks:
            if hasattr(self, 'undefined_bl') and bl in self.undefined_bl:
                gs.write_block(bl, self.undefined_bl[bl])
            else:
                do_warn('block name not found in undefined_bl:' + bl)

    def __write_general_block_v1(self, gs, bl, **kwargs):
        text2write = getattr(self, getattr(self._gr_block_names, bl))
        gs.write_block(bl, text2write)

    def __write_TITLE_v1(self, gs, bl, file_type = None, form = None, **kwargs):
        temp_lines = []
        first_line = []
        temp_cont = self.get_container(DescriptionPart, flag_class=True)
        for i in temp_cont:
            if i.file_type is not None and file_type is not None and i.file_type!=file_type:
                continue
            if i.form is not None and form is not None and i.form!=form:
                continue
            temp_lines.extend(first_line)
            first_line = ('', '')
            temp_lines.extend(i.lines)
        gs.write_block(bl, temp_lines)

    def __write_gromos_format_v1(self, gs, *blocks, **kwargs):
        check_unknown_bl = kwargs.get('check_unknown_bl', True)
        temp_s = None
        if kwargs.get('get_str', False):
            temp_s = ''
        for bl in blocks:
            if hasattr(self, '_write_' + bl):
                temp_f = getattr(self, getattr(self, '_write_' + bl))
            elif hasattr(self._gr_block_names, bl):
                temp_f = getattr(self, self.__write_general_block_v)
            elif check_unknown_bl:
                temp_f = self.__write_unknown_blocks
            temp = temp_f(gs, bl, **kwargs)
            if kwargs.get('get_str', False):
                temp_s += str(temp)
        return temp_s

    def write_gromos_format(self, gs, *blocks, **kwargs):
        temp_f = getattr(self, self.__write_gromos_format_v)
        return temp_f(gs, *blocks, **kwargs)

_GromosWriter_defs = {}
_GromosWriter_defs['__write_general_block_v'] = '__write_general_block_v1'
_GromosWriter_defs['__write_gromos_format_v'] = '__write_gromos_format_v1'
_GromosWriter_defs['_write_TITLE'] = '__write_TITLE_v1'

GromosWriter._add_defaults(_GromosWriter_defs, flag_set=True)

# FF (IFP) parsing & writing functions
block_interaction_type_map = dict()
block_interaction_type_map['BONDSTRETCHTYPECODE'] = BondType
block_interaction_type_map['BONDANGLEBENDTYPECODE'] = AngleType
block_interaction_type_map['IMPDIHEDRALTYPECODE'] = ImproperType
block_interaction_type_map['TORSDIHEDRALTYPECODE'] = DihedralType

top_block_interaction_type_map = dict()
for temp_bl in block_interaction_type_map:
    top_block_interaction_type_map[temp_bl[:-4]] = block_interaction_type_map[temp_bl]

class IFPBlocksParser(GromosParser):
    """contains functions that read of different IFP-related blocks"""
    def parse_ifp(self, parse_from, parse_from_file = True, **kwargs):
        """parses an ifp file"""
        self.parse_gro(parse_from, parse_from_file, **kwargs)
        self.find_imp_pairs()

    def __just_for_readability(self):
        self.generate_vdw = None # from FF
        self.add_vdw = None # from FF
        self.get_intDB = None  # from FF
        self.__read_interaction_type_code_v = None # from defaults defined below
        self.find_imp_pairs = None
        self._set_gr_interaction_type_container = None

    def __MASSATOMTYPECODE_v1(self, gs, *args, **kwargs):
        n = int(next(gs.block_split_fnc))
        n_last = int(next(gs.block_split_fnc))
        for i in range(n):
            m_type_id = next(gs.block_split_fnc)
            temp_masstype = MassType(m_type_id, next(gs.block_split_fnc), next(gs.block_split_fnc))
            self.get_intDB().add2container(temp_masstype, create = True, **kwargs)
        gs.end_block_split()

    def __SINGLEATOMLJPAIR_v1(self, gs, *args, **kwargs):
        read_a_type = []
        n_at = int(next(gs.block_split_fnc))
        for i in range(n_at):
            atom_id = next(gs.block_split_fnc)
            read_a_type.append(atom_id)
            atom_name = next(gs.block_split_fnc)
            temp_c612= []
            for j in range(6):
                temp_c612.append(float(next(gs.block_split_fnc)))
            temp_rules = []
            for j in range(n_at):
                temp_rules.append(int(next(gs.block_split_fnc)))
            temp_atom_type = AtomType(atom_id, atom_name, temp_c612, temp_rules, format_type='gr')
            self.get_intDB().add2container(temp_atom_type, create = True)
            if 'DUM' in atom_name:
                self.get_intDB().DUM_type = temp_atom_type
        self.generate_vdw(read_a_type)
        gs.end_block_split()

    def __MIXEDATOMLJPAIR_v1(self, gs, *args, **kwargs):
        flag = True
        while flag:
            try:
                at1 = next(gs.block_split_fnc)
            except StopIteration:
                return
            at2 = next(gs.block_split_fnc)
            c6 = next(gs.block_split_fnc)
            c12 = next(gs.block_split_fnc)
            c6_14 = next(gs.block_split_fnc)
            c12_14 = next(gs.block_split_fnc)
            vdw_normal = vdWType(params=[c6, c12])
            vdw_pair = PairType(params=[c6_14, c12_14])
            self.add_vdw(at1, at2, vdw_normal, vdw_pair, replace=True)

    def _read_block_interaction_type_code(self, gs, bl_name, **kwargs):
        temp_f = getattr(self, self.__read_interaction_type_code_v)
        interaction_type = block_interaction_type_map[bl_name]
        temp_f(interaction_type, gs, flag_code = True, **kwargs)

    def _read_block_interaction_type(self, gs, bl_name, **kwargs): # for topology parser
        temp_f = getattr(self, self.__read_interaction_type_code_v)
        interaction_type = top_block_interaction_type_map[bl_name]
        temp_f(interaction_type, gs, flag_code=False, **kwargs)

    def _read_interaction_type_code(self, bonded_type, in_f, flag_code=True, *args, **kwargs):
        temp_f = getattr(self, self.__read_interaction_type_code_v)
        temp_f(bonded_type, in_f, flag_code, *args, **kwargs)

    def __read_interaction_type_code_v1(self, bonded_type, in_f, flag_code=True, *args, **kwargs):
        n_b = int(next(in_f.block_split_fnc))
        if flag_code:
            _ = int(next(in_f.block_split_fnc))
        for i in range(n_b):
            if flag_code:
                b_code = next(in_f.block_split_fnc)
            else:
                b_code = str(i + 1)
            temp_b = bonded_type(b_code)
            self._set_gr_interaction_type_container(temp_b)
            for j in range(len(temp_b.flag)):
                temp_b.add_params(next(in_f.block_split_fnc))
            self.get_intDB().add2container(temp_b, create = True)
        in_f.end_block_split()

_GromosIFPBlocksParser_defs = {}
_GromosIFPBlocksParser_defs['_MASSATOMTYPECODE_parser'] = '__MASSATOMTYPECODE_v1'
_GromosIFPBlocksParser_defs['_SINGLEATOMLJPAIR_parser'] = '__SINGLEATOMLJPAIR_v1'
_GromosIFPBlocksParser_defs['_MIXEDATOMLJPAIR_parser'] = '__MIXEDATOMLJPAIR_v1'
_GromosIFPBlocksParser_defs['__read_interaction_type_code_v'] = '__read_interaction_type_code_v1'
for temp_bl in block_interaction_type_map:
    _GromosIFPBlocksParser_defs['_' + temp_bl + '_parser'] = '_read_block_interaction_type_code'
IFPBlocksParser._add_defaults(_GromosIFPBlocksParser_defs, flag_set=1)


class IFPBlocksWriter(GromosWriter):
    # writing part
    def write_ifp(self, f_path = None, **kwargs):
        gs = self._get_grs_from_fpaht(f_path)
        int_blocks2write = ('TITLE', 'FORCEFIELD', 'MAKETOPVERSION', 'MASSATOMTYPECODE', 'BONDSTRETCHTYPECODE',
                            'BONDANGLEBENDTYPECODE', 'IMPDIHEDRALTYPECODE', 'TORSDIHEDRALTYPECODE', 'SINGLEATOMLJPAIR',
                            'MIXEDATOMLJPAIR')
        self.write_gromos_format(gs, *int_blocks2write)
        if f_path and kwargs.get('flag_close', True):
            gs.f.close()

    def __just_for_readability(self):
        self.__form_a_type = tuple() # set in defaults
        self.__form_m_type = tuple() # set in defaults
        self.__form_interactions = tuple()  # set in defaults
        self.__write_interaction_type_v = None # set in defaults
        self.get_intDB = None

    def __find_last_m_code(self):
        m = 0
        for i in self.get_intDB().m_type:
            temp = int(i)
            if temp > m:
                m = temp
        return m

    def __write_MASSATOMTYPECODE_v1(self, gf, bl):
        m_type_cont = self.get_intDB().get_container('m_type', create = True)
        s = self.__form_m_type[0].format(len(m_type_cont), self.__find_last_m_code())
        for m_id, m_type in m_type_cont.items():
            s += self.__form_m_type[1].format(m_type.id, m_type.m, m_type.at_id)
        gf.write_block(bl, s)

    def __first_line_interaction(self, n, last_n=None):
        if not last_n:
            last_n = n
        return self.__form_interactions[0].format(n, last_n)

    def __write_interaction_type_v1(self, bonded_type, gf, bl, *args, **kwargs):
        temp_dict = self.get_intDB(bonded_type(None), cont_pref = self._gr_cont_pref, create=0, allow_not_found=1)
        if not temp_dict:
            gf.write_block(bl, self.__first_line_interaction(0))
            return
        s = self.__first_line_interaction(len(temp_dict))
        for int_code in temp_dict:
            s += self.__form_interactions[1].format(int_code)
            for i in range(len(temp_dict[int_code].p)):
                p = temp_dict[int_code].p[i]
                if temp_dict[int_code].p_type[i] is int:
                    s += self.__form_interactions[2].format(p)
                else:
                    s += self.__form_interactions[3].format(p)
            s += '\n'
        gf.write_block(bl, s)

    def __write_interaction_type(self, gf, bl, *args, **kwargs):
        temp_f = getattr(self, self.__write_interaction_type_v)
        interaction_type = block_interaction_type_map[bl]
        return temp_f(interaction_type, gf, bl, *args, **kwargs)

    def __write_SINGLEATOMLJPAIR_v1(self, gf, bl, *args, **kwargs):
        n_rulesline = self.__form_a_type[0]
        if not hasattr(self.get_intDB(), 'a_type'):return
        s = self.__form_a_type[1].format(len(self.get_intDB().a_type))
        for at_id in self.get_intDB().a_type:
            at = self.get_intDB().a_type[at_id]
            temp_s = self.__form_a_type[2] * 2 + self.__form_a_type[3] * 4 + '\n' + self.__form_a_type[3] * 2 + '\n'
            s += temp_s.format(get_id_string(at, flag_id='gr_id'), at.name, *at.c612)
            temp_s = ''
            for i in range(len(at.rules) // n_rulesline):
                temp_s += self.__form_a_type[4] * n_rulesline + '\n'
            temp_s += self.__form_a_type[4] * (len(at.rules) % n_rulesline) + '\n'
            s += temp_s.format(*at.rules)
        gf.write_block(bl, s)

    def __write_MIXEDATOMLJPAIR_v1(self, gf, bl):
        if not hasattr(self.get_intDB(), 'a_type'):return
        temp_ff = self.__class__()
        temp_ff.get_intDB().a_type = self.get_intDB().a_type
        temp_ff.generate_vdw(format_type = 'gr')
        s = ''
        for at_p in temp_ff.vdw:
            p1 = temp_ff.vdw[at_p][0].p + temp_ff.vdw[at_p][1].p
            if self.vdw[at_p][0] is not None and self.vdw[at_p][1] is not None:
                p2 = self.vdw[at_p][0].p + self.vdw[at_p][1].p
                flag_vdw_none = False
            else:
                flag_vdw_none = True
            if flag_vdw_none or not check_if_eq(p1, p2):
                s += (self.__form_a_type[2] * 2 + self.__form_a_type[3] * 4 + '\n').format(
                    get_id_string(at_p[0], flag_id='gr_id'), get_id_string(at_p[1], flag_id='gr_id'), *p2)
        gf.write_block(bl, s)

_GromosIFPBlocksWriter_defs = {}
_GromosIFPBlocksWriter_defs['_write_MASSATOMTYPECODE'] = '__write_MASSATOMTYPECODE_v1'
_GromosIFPBlocksWriter_defs['__write_interaction_type_v'] = '__write_interaction_type_v1'
for temp_bl in block_interaction_type_map:
    _GromosIFPBlocksWriter_defs['_write_' + temp_bl] = '__write_interaction_type'

_GromosIFPBlocksWriter_defs['_write_SINGLEATOMLJPAIR'] = '__write_SINGLEATOMLJPAIR_v1'
_GromosIFPBlocksWriter_defs['_write_MIXEDATOMLJPAIR'] = '__write_MIXEDATOMLJPAIR_v1'

IFPBlocksWriter._add_defaults(_GromosIFPBlocksWriter_defs, flag_set=1)

_GromosIFPBlocksWriter_cdefs = {}
_GromosIFPBlocksWriter_cdefs['form_m_type'] = ('{:>9d}{:>7d}\n', '{:<8}{:<7.4f} {:}\n')
_GromosIFPBlocksWriter_cdefs['form_interactions'] = ('{:9d}{:8d}\n', '{:>4}', '{:16d}', '{:16.7e}')
#first line, int_code, param(int), param(float)
_GromosIFPBlocksWriter_cdefs['form_a_type'] = (20, '{:10d}\n', '{:>8}', '{:14.5e}', '{:4d}')
# number_of_rules_inline, number of at_types, id and name, vdw params, rules
IFPBlocksWriter._add_class_defaults(_GromosIFPBlocksWriter_cdefs, flag_set=1)

# MTB and topology parsing and writing

# helper function for BB and TOP
def _parse_single_interaction(self, parse_from, interaction_type, ff, **kwargs):
    temp_interaction = Interaction(interaction_type)
    for _ in range(temp_interaction.na):
        temp_at = self.get_item(next(parse_from.block_split_fnc), self.Atom, create=True)
        temp_at.gr_id = temp_at.id
        temp_interaction.add_atom(temp_at)
    params_id = next(parse_from.block_split_fnc)
    temp_int_type = None
    if ff:
        temp_int_type = ff.get_intDB().get_item(params_id, interaction_type, allow_not_found=True,
                                                cont_pref = ff._gr_cont_pref)
        if temp_int_type:
            temp_interaction.add_state(temp_int_type)
    if not temp_int_type:
        temp_interaction.add_state(params_id)
    self.add2container(temp_interaction, db_type=list, create=True)


class BBParsing(GromosParser):
    def __just_for_readability(self):
        self.ChargeGroup = None
        self.sort_atoms = None
        self.Atom = None

    def _parse_bb_call(self, gs, ff = None, **kwargs):
        self.id = next(gs.block_split_fnc)
        n_at = int(next(gs.block_split_fnc))
        self.n_at = n_at
        n_pre_at = int(next(gs.block_split_fnc))
        self.n_pre_at = n_pre_at
        # preatoms
        for i in range(n_pre_at):
            self.__read_pre_atom(gs)
        # atoms
        cg = self.get_item(None, self.ChargeGroup, create=True, create_container=True, db_type=list)
        for i in range(n_at - n_pre_at):
            cg = self.__read_atom(gs, cg, ff)
        # postatoms
        for i in range(n_pre_at):
            cg = self.__read_atom(gs, cg, ff)
        self.cg.pop()
        self.__parse_interaction(gs, BondType, ff)
        self.__parse_interaction(gs, AngleType, ff)
        self.__parse_interaction(gs, ImproperType, ff)
        self.__parse_interaction(gs, DihedralType, ff)
        self.lj_exceptions = list(gs.block_lines())
        self.sort_atoms()

    def __MTBUILDBLSOLUTE_v1(self, gs, ff = None, **kwargs):
        self._parse_bb_call(gs, ff = ff, **kwargs)

    __parse_single_interaction = _parse_single_interaction

    def __parse_interaction(self, gs, interaction, ff, **kwargs):
        n_b = int(next(gs.block_split_fnc))
        for i in range(n_b):
            self.__parse_single_interaction(gs, interaction, ff, **kwargs)

    def __read_pre_atom(self, gs):
        at_id = next(gs.block_split_fnc)
        temp_at = self.get_item(at_id, self.Atom, create=True, create_container = True)
        temp_at.gr_id = at_id
        n_ex = int(next(gs.block_split_fnc))
        for i in range(n_ex):
            temp_at_n = next(gs.block_split_fnc)
            temp_at_e = self.get_item(temp_at_n, self.Atom, create=True)
            temp_at_e.gr_id = temp_at_n
            temp_at.add_excl(temp_at_e)

    def __read_atom(self, gs, cg, ff=None, **kwargs):
        at_id = next(gs.block_split_fnc)
        temp_at = self.get_item(at_id, self.Atom, create=True, create_container = True)
        temp_at.gr_id = at_id
        temp_at.name = next(gs.block_split_fnc)
        temp_at.a_type = next(gs.block_split_fnc)
        temp_at.m_type = next(gs.block_split_fnc)
        if ff:
            temp_true_type = ff.get_intDB().get_item(temp_at.a_type, AtomType, allow_not_found=True)
            if temp_true_type:
                temp_at.a_type = temp_true_type
            temp_true_type = ff.get_intDB().get_item(temp_at.m_type, MassType, allow_not_found=True)
            if temp_true_type:
                temp_at.m_type = temp_true_type
        temp_at.p_ch = float(next(gs.block_split_fnc))
        cg.add_atom(temp_at)
        mk_cg = int(next(gs.block_split_fnc))
        temp_at.mk_cg = mk_cg
        if mk_cg:
            #cg = self.ChargeGroup()
            cg = self.get_item(None, self.ChargeGroup, create=True, create_container=True, db_type=list)
        flag_excl = not gs.report_nl_block_split_fnc
        if flag_excl:
            n_excl = int(next(gs.block_split_fnc))
            for i in range(n_excl):
                temp_at_n = next(gs.block_split_fnc)
                temp_at_e = self.get_item(temp_at_n, self.Atom, create=True)
                temp_at_e.gr_id = temp_at_n
                temp_at.add_excl(temp_at_e)
        temp_at.flag_bb = True
        temp_at.flag_excl = flag_excl
        return cg

_GromosBBParsing_defs = {}
_GromosBBParsing_defs['_MTBUILDBLSOLUTE_parser'] = '__MTBUILDBLSOLUTE_v1'
BBParsing._add_defaults(_GromosBBParsing_defs, flag_set=1)


class grAtomWriting:
    gr_id_format = '{:>5}'
    n_atoms_format = '{:4d}'
    new_line_format = '\n' + ' ' * 47

    @property
    def write_gr_id(self):
        if hasattr(self, 'gr_id'):
            gr_id = getattr(self, 'gr_id')
        else:
            gr_id = getattr(self, 'id')
        return self.gr_id_format.format(gr_id)

    def _write_gr_EP(self, e_p_l):
        n = len(e_p_l)
        s = self.n_atoms_format.format(n)
        for i in range(n):
            if i % 6 == 0 and i != 0:
                s += self.new_line_format
            e_at = e_p_l[i]
            s += e_at.write_gr_id
        return s + '\n'


class grBBAtomWriting(grAtomWriting):
    gr_id_format = '{:>5}'
    n_atoms_format = '{:4d}'
    new_line_format = '\n' + ' ' * 39

    def __just_for_readability(self):
        self.flag_excl = None
        self.e_l = None

    @property
    def write_gr_excl(self):
        if not self.flag_excl:
            return '\n'
        """
        n = len(self.e_l)
        s = '{:4d}'.format(n)
        for i, e_at in enumerate(self.e_l):
            if i % 6 == 0 and i != 0:
                s += '\n' + ' ' * 39
            s += e_at.write_gr_id
        return s + '\n'
        """
        return self._write_gr_EP(self.e_l)

    @property
    def write_gr_prop(self):
        temp_at = get_id_string(self.a_type, flag_id='gr_id')
        temp_mt = get_id_string(self.m_type)
        s = ' {:<4}{:>5}{:>5}{:11.5f}{:>4}'.format(self.name, temp_at, temp_mt, self.p_ch, self.mk_cg)
        return s


class grTOPAtomWriting(grAtomWriting):
    gr_id_format = '{:>6}'
    n_atoms_format = '{:6d}'
    new_line_format = '\n' + ' ' * 47

    def __just_for_readability(self):
        self.id = None
        self.gr_id = None

    @staticmethod
    def __write_gr_EP(e_p_l):
        n = len(e_p_l)
        s = '{:6d}'.format(n)
        for i in range(n):
            if i % 6 == 0 and i != 0:
                s += '\n' + ' ' * 47
            e_at = e_p_l[i]
            s += e_at.write_gr_id
        return s + '\n'

    def write_gr_excl(self, e_l):
        return self._write_gr_EP(e_l)

    def write_gr_pairs(self, p_l):
        return ' ' * 41 + self._write_gr_EP(p_l)

    @property
    def write_gr_prop(self, rev_prop = False):
        if hasattr(self, 'at_format'):
            at_format = self.at_format
        else:
            at_format = '{:>5}{:>5}{:>4}{:9.5f}{:9.5f}{:>3d}'
        temp_at = get_id_string(self.a_type, flag_id='gr_id')
        if not rev_prop:
            at_name = self.name
        else:
            at_name = getattr(self, 'rev_name', self.name)
        s = at_format.format(get_id_string(self.res, flag_id='gr_id'), at_name, temp_at, self.m, self.p_ch, self.mk_cg)
        return s

    @property
    def write_rev_prop(self):
        return self.write_gr_prop(rev_prop = True)


class BBWriting(GromosWriter):
    def __write_MTBUILDBLSOLUTE_v1(self, gf, bl, **kwargs):
        s = self.id + '\n'
        bb_atoms = self.__get_bb_atoms()
        s += '{:5d}{:5d}\n'.format(len(bb_atoms), self.n_pre_at)
        s += '# preceding exclusions\n#ATOM                               MAE MSAE\n'
        for i in range(-1 * self.n_pre_at + 1, 1):
            pre_at = self.__get_pre_atom(i)
            s += pre_at.write_gr_id + ' ' * 29
            s += pre_at.write_gr_excl
        s += '# atoms\n# ATOM ANM  IACM MASS        CGMICGM MAE MSAE\n'
        for at in bb_atoms:
            s += at.write_gr_id + at.write_gr_prop + at.write_gr_excl

        for temp_int_type in (BondType, AngleType, ImproperType, DihedralType):
            s += '# ' + temp_int_type.container2write + '\n'
            temp_db = self.get_container(Interaction(temp_int_type), flag_item=True, create = True)
            s += '{:5d}\n'.format(len(temp_db))
            for temp_int in temp_db:
                for at in temp_int.atoms:
                    s += at.write_gr_id
                s += '{:>5}\n'.format(get_id_string(temp_int.states[0]))
        s += '# LJ exceptions\n'
        for i in self.lj_exceptions:
            s += i
        if kwargs.get('get_str', False):
            return s
        else:
            gf.write_block(bl, s)

    def __get_pre_atom(self, at_id):
        return self.atoms[str(at_id)]

    def __get_bb_atoms(self):
        bb_at = []
        for at in self.atoms:
            if self.atoms[at].flag_bb:
                bb_at.append(self.atoms[at])
        return bb_at

_GromosBBWriting_defs = {}
_GromosBBWriting_defs['_write_MTBUILDBLSOLUTE'] = '__write_MTBUILDBLSOLUTE_v1'
BBWriting._add_defaults(_GromosBBWriting_defs, flag_set=1)

class MTBBlocksParser(GromosParser):

    def __MTBUILDBLSOLUTE_v1(self, gs, **kwargs):
        bb = self.BuildingBlock()
        bb._parse_bb_call(gs, ff=self, **kwargs)
        self.add2container(bb, create=True)

    def parse_mtb(self, parse_from, parse_from_file = True, **kwargs):
        """parses a mtb file and set appropriate parameters within the instance of MTB"""
        self.parse_gro(parse_from, parse_from_file, **kwargs)



_GromosMTBBlocksParser_defs = {}# same as _GromosBBParsing_defs
_GromosMTBBlocksParser_defs['_MTBUILDBLSOLUTE_parser'] = '__MTBUILDBLSOLUTE_v1'
MTBBlocksParser._add_defaults(_GromosMTBBlocksParser_defs, flag_set=1)

### Topology
top_block_interaction_map = dict()
top_block_interaction_map['BONDH'] = BondType
top_block_interaction_map['BOND'] = BondType
top_block_interaction_map['BONDANGLEH'] = AngleType
top_block_interaction_map['BONDANGLE'] = AngleType
top_block_interaction_map['IMPDIHEDRALH'] = ImproperType
top_block_interaction_map['IMPDIHEDRAL'] = ImproperType
top_block_interaction_map['DIHEDRALH'] = DihedralType
top_block_interaction_map['DIHEDRAL'] = DihedralType
top_block_interaction_map['CONSTRAINT'] = ConstraintType

ptp_block_interaction_map = dict()
ptp_block_interaction_map['PERTBONDSTRETCH'] = BondType
ptp_block_interaction_map['PERTBONDANGLE'] = AngleType
ptp_block_interaction_map['PERTIMPROPERDIH'] = ImproperType
ptp_block_interaction_map['PERTPROPERDIH'] = DihedralType


class topBlocksParser(GromosParser):

    def __just_for_readability(self):
        self.sort_atoms = None
        self.get_pairs_from_p_l = None
        self.add_a_type = None
        self.add_vdw = None
        self.get_intDB = None
        self.get_a_type = None
        self.add_residue = None
        self.Atom = None
        self.Residue = None
        self.ChargeGroup = None

    def parse_top_gr(self, parse_from, parse_from_file = True, **kwargs):
        """parses a top file (GROMOS format)"""
        self.parse_gro(parse_from, parse_from_file, **kwargs)
        self.cg.pop()
        self.sort_atoms()
        #self.get_pairs_from_p_l()

    def __ATOMTYPENAME_v1(self, in_f, *args, **kwargs):
        n_entries = int(next(in_f.block_split_fnc))
        for i in range(n_entries):
            atom_name = next(in_f.block_split_fnc)
            a_type_id = str(i + 1)
            self.add_a_type(a_type_id, atom_name, format_type='gr')
            if 'DUM' in atom_name:
                self.get_intDB().DUM_type = self.get_a_type(a_type_id)
        in_f.end_block_split()

    def __RESNAME_v1(self, in_f, *args, **kwargs):
        n_entries = int(next(in_f.block_split_fnc))
        for i in range(n_entries):
            self.add_residue(res_name=next(in_f.block_split_fnc))
        in_f.end_block_split()

    def __SOLUTEATOM_v1(self, in_f, *args, **kwargs):
        n_entries = int(next(in_f.block_split_fnc))
        cg = self.get_item(None, self.ChargeGroup, create=True, create_container=True, db_type = list)
        for _ in range(n_entries):
            at = self.get_item(next(in_f.block_split_fnc), self.Atom, create=True, create_container=True)
            at.gr_id = at.id
            res = self.get_item(next(in_f.block_split_fnc), self.Residue, create_container=True)
            res.gr_id = res.id
            res.add_atom(at)
            #at.res = res
            at.name = next(in_f.block_split_fnc)
            at.a_type = self.get_intDB().a_type[next(in_f.block_split_fnc)]
            at.m = float(next(in_f.block_split_fnc))
            at.p_ch = float(next(in_f.block_split_fnc))
            at._generate_self_state()
            cg.add_atom(at)
            mk_cg = int(next(in_f.block_split_fnc))
            at.mk_cg = mk_cg
            if mk_cg:
                cg = self.get_item(None, self.ChargeGroup, create=True, create_container=True, db_type=list)
            n_excl = int(next(in_f.block_split_fnc))
            for _ in range(n_excl):
                temp_at = self.get_item(next(in_f.block_split_fnc), self.Atom, create=True, create_container=True)
                #at.add_excl(temp_at)
                self.add_atom_pair2EP_l(at, temp_at)
                """
                temp_interaction = Interaction(ExclusionType, atoms=(at, temp_at))
                temp_interaction.add_state(fnc_type = 'gr_fnc', params = (True,))
                #self.add2container(temp_interaction, create=True, db_type=list)
                """
                # exclusion pair type
                temp_interaction = Interaction(ExclusionPairType, atoms=(at, temp_at))
                temp_interaction.add_state(fnc_type = 'gr_fnc', params = (True,))
                self.add2container(temp_interaction, create=True, item_id=frozenset((at, temp_at)), replace = -1)
            n_pairs = int(next(in_f.block_split_fnc))
            for _ in range(n_pairs):
                temp_at = self.get_item(next(in_f.block_split_fnc), self.Atom, create=True, create_container=True)
                #at.add_pair(temp_at)
                self.add_atom_pair2EP_l(at, temp_at)
                temp_interaction = Interaction(ExclusionPairType, atoms=(at, temp_at))
                temp_interaction.add_state(fnc_type = 'gr_fnc', params = (2,))
                self.add2container(temp_interaction, create=True, item_id=frozenset((at, temp_at)), replace = -1)
                #temp_atoms = (at, temp_at)
                #temp_state = self.find_interaction_type(temp_atoms, PairType, 'gr_fnc', **kwargs)
        in_f.end_block_split()

    __parse_single_interaction = _parse_single_interaction

    def _parse_interaction(self, gs, bl_name, **kwargs):
        interaction_type = top_block_interaction_map[bl_name]
        self.__parse_interaction(gs, interaction_type, self, **kwargs)

#    def parse_interaction(self, parse_from, interaction, ff=None, **kwargs):
#        self.__parse_interaction(parse_from, interaction, ff, **kwargs)

    def __parse_interaction(self, parse_from, interaction, ff, **kwargs):
        n_b = int(next(parse_from.block_split_fnc))
        for i in range(n_b):
            self.__parse_single_interaction(parse_from, interaction, ff, **kwargs)
        parse_from.end_block_split()

    def __LJPARAMETERS_v1(self, parse_from, *args, **kwargs):
        n_b = int(next(parse_from.block_split_fnc))
        for i in range(n_b):
            at1, at2 = [next(parse_from.block_split_fnc) for _ in range(2)]
            vdw_normal = vdWType(params = [next(parse_from.block_split_fnc) for _ in range(2)][::-1])
            vdw_pair = PairType(params = [next(parse_from.block_split_fnc) for _ in range(2)][::-1])
            vdw_normal.convert2se()
            vdw_pair.convert2se()
            self.add_vdw(at1, at2, vdw_normal, vdw_pair)
        parse_from.end_block_split()

    def __read_groups(self, parse_from, cont_name):
        cont_name = 'db_' + cont_name
        n = int(next(parse_from.block_split_fnc))
        if not hasattr(self, cont_name):
            setattr(self, cont_name, list())
        cont = getattr(self, cont_name)
        for _ in range(n):
            cont.append(next(parse_from.block_split_fnc))
        parse_from.end_block_split()

    def __SOLUTEMOLECULES_v1(self, parse_from, bl_name):
        self.__read_groups(parse_from, bl_name)

    def __TEMPERATUREGROUPS_v1(self, parse_from, bl_name):
        self.__read_groups(parse_from, bl_name)

    def __PRESSUREGROUPS_v1(self, parse_from, bl_name):
        self.__read_groups(parse_from, bl_name)


    def parse_dsr(self, parse_from, parse_from_file = True, **kwargs):
        """parses a top file (GROMOS format)"""
        self.parse_gro(parse_from, parse_from_file, **kwargs)

    def __DISTANCERESSPEC_v1(self, parse_from, bl_name):
        dish, disc = float(next(parse_from.block_split_fnc)), float(next(parse_from.block_split_fnc))
        if hasattr(self, 'gr_distres_dis_const'):
            np.testing.assert_almost_equal(self.gr_distres_dis_const, [dish, disc])
        else:
            self.gr_distres_dis_const = [dish, disc]
        for l in parse_from.block_lines():
            temp_dr = Interaction(self.Distance_r)
            temp = split_gen(l)
            temp_atoms = []
            temp_dr.gr_extra_atoms = [[], []]
            temp_dr.gr_atom_types = []
            for i in range(2):
                temp_atoms.append(self.atoms[next(temp)])
                for j in range(3):
                    temp_at = next(temp)
                    if temp_at == '0':
                        temp_dr.gr_extra_atoms[i].append(None)
                    else:
                        temp_dr.gr_extra_atoms[i].append(self.atoms[temp_at])
                temp_dr.gr_atom_types.append(next(temp))
            temp_dr.add_atom(*temp_atoms)
            self.add2container(temp_dr, create=True, db_type=list)
            temp_dr.add_state()
            temp_dr.states[-1].add_params(temp)

_topBlocksParser_defs = {}
_topBlocksParser_defs['_ATOMTYPENAME_parser'] = '__ATOMTYPENAME_v1'
_topBlocksParser_defs['_RESNAME_parser'] = '__RESNAME_v1'
_topBlocksParser_defs['_SOLUTEATOM_parser'] = '__SOLUTEATOM_v1'
for temp_bl in top_block_interaction_type_map:
    _topBlocksParser_defs['_' + temp_bl + '_parser'] = '_read_block_interaction_type'
for temp_bl in top_block_interaction_map:
    _topBlocksParser_defs['_' + temp_bl + '_parser'] = '_parse_interaction'
_topBlocksParser_defs['_LJPARAMETERS_parser'] = '__LJPARAMETERS_v1'
_topBlocksParser_defs['_SOLUTEMOLECULES_parser'] = '__SOLUTEMOLECULES_v1'
_topBlocksParser_defs['_TEMPERATUREGROUPS_parser'] = '__TEMPERATUREGROUPS_v1'
_topBlocksParser_defs['_PRESSUREGROUPS_parser'] = '__PRESSUREGROUPS_v1'
_topBlocksParser_defs['_DISTANCERESSPEC_parser'] = '__DISTANCERESSPEC_v1'
topBlocksParser._add_defaults(_topBlocksParser_defs, flag_set=True)


class topBlocksWriter(GromosWriter):

    def __just_for_readability(self):
        self.get_intDB = None # from FF
        self.residues = OrderedDict()
        self.atoms = OrderedDict()
        self._gr_cont_pref = None # from FF
        self.__write_interaction_atoms_v = None # defaults
        self.__write_interaction_type_v = None # defaults
        self._gr_molecules_string = '' # from Topology (GROMOS format)
        self.get_molecules = None # from Topology (GROMOS format)

    def write_top_gr(self, f_path = None, **kwargs):
        self._get_gr_molecule_string()
        gs = self._get_grs_from_fpaht(f_path)
        self.write_gromos_format(gs, 'TITLE', 'PHYSICALCONSTANTS', 'TOPVERSION',)
        if kwargs.get('constraints', False):
            int_blocks2write = ('ATOMTYPENAME', 'RESNAME', 'SOLUTEATOM', 'BONDSTRETCHTYPE', 'BONDH', 'BOND',
                                'CONSTRAINT', 'BONDANGLEBENDTYPE', 'BONDANGLEH', 'BONDANGLE', 'IMPDIHEDRALTYPE',
                                'IMPDIHEDRALH', 'IMPDIHEDRAL', 'TORSDIHEDRALTYPE', 'DIHEDRALH', 'DIHEDRAL',
                                'LJPARAMETERS', 'SOLUTEMOLECULES', 'TEMPERATUREGROUPS', 'PRESSUREGROUPS',
                                'LJEXCEPTIONS', 'SOLVENTATOM', 'SOLVENTCONSTR')
        else:
            int_blocks2write = ('ATOMTYPENAME', 'RESNAME', 'SOLUTEATOM', 'BONDSTRETCHTYPE', 'BONDH', 'BOND',
                                'BONDANGLEBENDTYPE', 'BONDANGLEH', 'BONDANGLE', 'IMPDIHEDRALTYPE',
                                'IMPDIHEDRALH', 'IMPDIHEDRAL', 'TORSDIHEDRALTYPE', 'DIHEDRALH', 'DIHEDRAL',
                                'LJPARAMETERS', 'SOLUTEMOLECULES', 'TEMPERATUREGROUPS', 'PRESSUREGROUPS',
                                'LJEXCEPTIONS', 'SOLVENTATOM', 'SOLVENTCONSTR')
        self.write_gromos_format(gs, *int_blocks2write, **kwargs)
        if f_path and kwargs.get('flag_close', True):
            gs.f.close()
# missing CROSSDIHEDRALH, CROSSDIHEDRAL

    def write_ptp(self, f_path, **kwargs):
        self._get_gr_molecule_string()
        temp_f = GromosFile(f_path, write=True)
        self.write_gromos_format(temp_f, 'TITLE')
        int_blocks2write = ('PERTATOMPARAM', 'PERTATOMPAIR', 'PERTBONDSTRETCH',
                            'PERTBONDANGLE', 'PERTIMPROPERDIH', 'PERTPROPERDIH')
        top_state, other_state = self._get_top_other_state(**kwargs)
        top_other_state = dict(top_state=top_state, other_state=other_state)
        write_kwargs = dict(kwargs)
        write_kwargs.update(top_other_state)
        self.write_gromos_format(temp_f, *int_blocks2write, **write_kwargs)

    def write_rev_ptp_top_names(self, f_path, **kwargs):
        self.write_top_gr(f_path, rev_top=True, **kwargs)

    def __write_ATOMTYPENAME_v1(self, gf, bl, *args, **kwargs):
        s = str(len(self.get_intDB().a_type)) + '\n'
        for at in self.get_intDB().a_type:
            s += self.get_intDB().a_type[at].name + '\n'
        gf.write_block(bl, s)

    def __write_RESNAME_v1(self, gf, bl, *args, **kwargs):
        s = str(len(self.residues)) + '\n'
        for res in self.residues:
            if kwargs.get('rev_top', False):
                rev_res_name = getattr(self.residues[res], 'rev_name', False)
                if rev_res_name:
                    s += rev_res_name + '\n'
                else:
                    s += self.residues[res].name + '\n'
            else:
                s += self.residues[res].name + '\n'
        gf.write_block(bl, s)

    def __write_SOLUTEATOM_v1(self, gf, bl, *args, **kwargs):
        e_l, p_l = self._generate_atom_excl_pair_list(**kwargs)
        s = str(len(self.atoms)) + '\n'
        for at in self.get_atoms():
            s += at.write_gr_id
            if kwargs.get('rev_top', False):
                s += at.write_rev_prop
            else:
                s += at.write_gr_prop
            s += at.write_gr_excl(e_l[at])
            s += at.write_gr_pairs(p_l[at])
        gf.write_block(bl, s)

    def __write_interaction_type_v1(self, bonded_type, gf, bl, *args, **kwargs):
        temp_dict = self.get_intDB(bonded_type, cont_pref = self._gr_cont_pref, create = False, allow_not_found = True)
        if not temp_dict:
            temp_dict = {}
        s = str(len(temp_dict)) + '\n'
        for int_code in temp_dict:
            for i in range(len(temp_dict[int_code].p)):
                p = temp_dict[int_code].p[i]
                if temp_dict[int_code].p_type[i] is int:
                    s += '{:4d}'.format(p)
                else:
                    s += '{:16.7e}'.format(p)
            s += '\n'
        gf.write_block(bl, s)

    def _write_interaction_type(self, gf, bl, *args, **kwargs):
        temp_f = getattr(self, self.__write_interaction_type_v)
        interaction_type = top_block_interaction_type_map[bl]
        return temp_f(interaction_type, gf, bl, *args, **kwargs)

    def __write_interaction_atoms_v1(self, bonded_container, gf, bl, *args, **kwargs):
        if bonded_container is None:
            return
        top_state = kwargs.get('top_state', getattr(self, 'top_state', 0))
        s = str(len(bonded_container)) + '\n'
        for bond in bonded_container:
            for at in bond.atoms:
                s += ' ' + at.write_gr_id
            state4writing = self._get_state(bond, **kwargs)
            """
            if hasattr(bond, 'state'):
                state4writing = bond.state
            else:
                state4writing = self.get_state(bond.states, **kwargs)
            """
            """
            if state4writing is None:
                raise Exception('all interaction states None: ' + str(bond.int_type) + '\natoms: ' + str(bond.atoms))
            if top_state == 'any_non_none':
                flag_state_not_found = True
                for temp_state in bond.states:
                    if temp_state is not None:
                        s += '{:>5}\n'.format(get_id_string(temp_state))
                        flag_state_not_found = False
                        break
                if flag_state_not_found:
                    raise Exception('all interaction states None: ' + str(bond.int_type) + '\natoms: ' + str(bond.atoms))
            else:
                s += '{:>5}\n'.format(get_id_string(bond.states[top_state]))
            """
            s += '{:>5}\n'.format(get_id_string(state4writing))
        gf.write_block(bl, s)

    def _write_interaction_atoms(self, gf, bl, *args, **kwargs):
        temp_f = getattr(self, self.__write_interaction_atoms_v)
        interaction_type = top_block_interaction_map[bl]
        bonded_container = self.get_container(interaction_type, flag_class = 1, allow_not_found = 1)
        return temp_f(bonded_container, gf, bl, *args, **kwargs)

    def __write_empty_interaction_atoms(self, gf, bl, *args, **kwargs):
        gf.write_block(bl, '0')

    def __write_LJPARAMETERS_v1(self, gf, bl, *args, **kwargs):
        s = str(len(self.get_intDB().vdw)) + '\n'
        for at_p in self.get_intDB().vdw:
            vdw_p = self.get_intDB().vdw[at_p]
            s += '{:>5}{:>5}'.format(get_id_string(at_p[0], flag_id='gr_id'), get_id_string(at_p[1], flag_id='gr_id'))
            temp_str = '{:14.6e}' * 4 + '\n'
            if vdw_p[1] is None:
                pair_p = (0,0)
            else:
                pair_p =  vdw_p[1].c612
#            s += temp_str.format(vdw_p[0].p[1], vdw_p[0].p[0], vdw_p[1].p[1], vdw_p[1].p[0])
            s += temp_str.format(vdw_p[0].c612[1], vdw_p[0].c612[0], pair_p[1], pair_p[0])
        gf.write_block(bl, s)

    def __write_mol_string_v1(self, gf, bl, *args, **kwargs):
        gf.write_block(bl, self._gr_molecules_string)

    def _get_gr_molecule_string(self):
        mols = self.get_molecules(flag_mol_atom_sort = 1)
        s = str(len(mols)) + '\n'
        for i,m in enumerate(mols):
            if i % 10 == 0 and i != 0:
                s += '\n'
            s += ' {:>5}'.format(get_id_string(m[-1], flag_id='gr_id'))
        self._gr_molecules_string = s

# ptp writing
    def __write_PERTATOMPARAM_v1(self, gf, bl, *args, **kwargs):
        DUM_type = self.get_DUM_type
        ALJ = kwargs.get('ALJ', 1.0)
        ACRF = kwargs.get('ACRF', 1.0)
        #top_state, other_state = self._get_top_other_state(**kwargs)
        a_type_kwargs, m_kwargs, pch_kwargs = self._get_atom_kwargs(**kwargs)
        s = ''
        N_at_ptp = 0
        for at in self.get_atoms():
            a_type_ptp = self.get_ptp_states(at, **a_type_kwargs)
            m_ptp = self.get_ptp_states(at, **m_kwargs)
            p_ch_ptp = self.get_ptp_states(at, **pch_kwargs)
            a_type_ptp, m_ptp, p_ch_ptp, flag_ptp = self.get_atom_ptp_states(at, a_type_ptp, m_ptp, p_ch_ptp, **kwargs)
            if flag_ptp:
                N_at_ptp += 1
                s += at.write_gr_id
                at_prop = [get_id_string(at.res, flag_id='gr_id'), at.name]
                ptp_prop = []
                at_prop.append(get_id_string(a_type_ptp[0], flag_id='gr_id'))
                ptp_prop.append(get_id_string(a_type_ptp[1], flag_id='gr_id'))
                at_prop.append(m_ptp[0])
                ptp_prop.append(m_ptp[1])
                at_prop.append(p_ch_ptp[0])
                ptp_prop.append(p_ch_ptp[1])
                # get the text
                s += '{:>5}{:>6}{:>4}{:11.5f}{:11.5f}'.format(*at_prop)
                s += '{:>4}{:11.5f}{:11.5f}'.format(*ptp_prop)
                # get the softness
                if not kwargs.get('flag_soft_heawy', False) and \
                        (m_ptp[0] > 1.5 and m_ptp[0] > 1.5 and a_type_ptp[0] != DUM_type and a_type_ptp[0] != DUM_type):
                    s += '{:9d}{:9d}\n'.format(0, 0)
                else:
                    s += '{:9.2f}{:9.2f}\n'.format(ALJ, ACRF)
        s = str(N_at_ptp) + '\n' + s
        gf.write_block(bl, s)

    def __write_PERTATOMPAIR_v1(self, gf, bl, *args, **kwargs):
        N_excl_pair_ptp = 0
        s = ''
        for at_p in self.excl_pair:
            excl_pair_ptp = self.get_ptp_states(self.excl_pair[at_p], **kwargs)
            if excl_pair_ptp:
                for at in at_p:
                    s += at.write_gr_id
                for i in range(2):
                    s += '{:>6}'.format(excl_pair_ptp[i]._code_map[excl_pair_ptp[i].p[0]])
                s += '\n'
                N_excl_pair_ptp += 1
        s = str(N_excl_pair_ptp) + '\n' + s
        gf.write_block(bl, s)

    def __write_pert_interaction_v1(self, bonded_container, gf, bl, *args, **kwargs):
        if bonded_container is None:
            return
        N_int_ptp = 0
        s = ''
        for temp_int in bonded_container:
            temp_int_ptp = self.get_ptp_states(temp_int)
            if temp_int_ptp:
                N_int_ptp += 1
                for at in temp_int.atoms:
                    s += ' ' + at.write_gr_id
                for i in range(2):
                    s += '{:>6}'.format(get_id_string(temp_int_ptp[i]))
                s += '\n'
        if N_int_ptp==0:
            return
        s = str(N_int_ptp) + '\n' + s
        return gf.write_block(bl, s)

    def _write_pert_interaction(self, gf, bl, *args, **kwargs):
        temp_f = getattr(self, self.__write_pert_interaction_v)
        interaction_type = ptp_block_interaction_map[bl]
        bonded_container = self.get_container(interaction_type, flag_class = 1, allow_not_found = 1)
        temp_f(bonded_container, gf, bl, *args, **kwargs)

    def __write_PERTBONDSTRETCH_v1(self, gf, bl, *args, **kwargs):
        self._write_pert_interaction(gf, bl, self.bonds, *args, **kwargs)

    def __write_PERTBONDANGLE_v1(self, gf, bl, *args, **kwargs):
        self._write_pert_interaction(gf, bl, self.angles, *args, **kwargs)

    def __write_PERTIMPROPERDIH_v1(self, gf, bl, *args, **kwargs):
        self._write_pert_interaction(gf, bl, self.impropers, *args, **kwargs)

    def __write_PERTPROPERDIH_v1(self, gf, bl, *args, **kwargs):
        self._write_pert_interaction(gf, bl, self.dihedrals, *args, **kwargs)


_topBlocksWriter_defs = {}
_topBlocksWriter_defs['_write_ATOMTYPENAME'] = '__write_ATOMTYPENAME_v1'
_topBlocksWriter_defs['_write_RESNAME'] = '__write_RESNAME_v1'
_topBlocksWriter_defs['_write_SOLUTEATOM'] = '__write_SOLUTEATOM_v1'
_topBlocksWriter_defs['__write_interaction_type_v'] = '__write_interaction_type_v1'

for temp_bl in top_block_interaction_type_map:
    _topBlocksWriter_defs['_write_' + temp_bl] = '_write_interaction_type'

_topBlocksWriter_defs['_write_SINGLEATOMLJPAIR'] = '__write_SINGLEATOMLJPAIR_v1'
_topBlocksWriter_defs['_write_MIXEDATOMLJPAIR'] = '__write_MIXEDATOMLJPAIR_v1'
_topBlocksWriter_defs['__write_interaction_atoms_v'] = '__write_interaction_atoms_v1'

for temp_bl in top_block_interaction_map:
    _topBlocksWriter_defs['_write_' + temp_bl] = '_write_interaction_atoms'

_topBlocksWriter_defs['_write_BONDH'] = '__write_empty_interaction_atoms'
_topBlocksWriter_defs['_write_BONDANGLEH'] = '__write_empty_interaction_atoms'
_topBlocksWriter_defs['_write_IMPDIHEDRALH'] = '__write_empty_interaction_atoms'
_topBlocksWriter_defs['_write_DIHEDRALH'] = '__write_empty_interaction_atoms'

_topBlocksWriter_defs['_write_LJPARAMETERS'] = '__write_LJPARAMETERS_v1'
_topBlocksWriter_defs['_write_SOLUTEMOLECULES'] = '__write_mol_string_v1'
_topBlocksWriter_defs['_write_TEMPERATUREGROUPS'] = '__write_mol_string_v1'
_topBlocksWriter_defs['_write_PRESSUREGROUPS'] = '__write_mol_string_v1'

_topBlocksWriter_defs['__write_pert_interaction_v'] = '__write_pert_interaction_v1'
_topBlocksWriter_defs['_write_PERTATOMPARAM'] = '__write_PERTATOMPARAM_v1'
_topBlocksWriter_defs['_write_PERTATOMPAIR'] = '__write_PERTATOMPAIR_v1'

for temp_bl in ptp_block_interaction_map:
    _topBlocksWriter_defs['_write_' +  temp_bl] = '_write_pert_interaction'

"""
_topBlocksWriter_defs['_write_PERTBONDSTRETCH'] = '__write_PERTBONDSTRETCH_v1'
_topBlocksWriter_defs['_write_PERTBONDANGLE'] = '__write_PERTBONDANGLE_v1'
_topBlocksWriter_defs['_write_PERTIMPROPERDIH'] = '__write_PERTIMPROPERDIH_v1'
_topBlocksWriter_defs['_write_PERTPROPERDIH'] = '__write_PERTPROPERDIH_v1'
"""

topBlocksWriter._add_defaults(_topBlocksWriter_defs, flag_set=1)


"""
class ptpBlocksWriter(GromosWriter): ####################################### something to change!!!

    def __write_PERTATOMPARAM_v1(self, gf, bl, *args, **kwargs):
        ALJ = kwargs.get('ALJ', 1.0)
        ACRF = kwargs.get('ACRF', 1.0)
        s = str(len(self.ptp.atoms)) + '\n'
        for at in self.ptp.atoms:
            at_ptp = self.ptp.atoms[at]
            temp_at = get_id_string(at.a_type)
            s += at.write_gr_id
            s += '{:>5}{:>6}{:>4}{:11.5f}{:11.5f}'.format(get_id_string(at.res, flag_id='gr_id'),
                                                          at.name, temp_at, at.m, at.p_ch)
            s += '{:>4}{:11.5f}{:11.5f}'.format(get_id_string(at_ptp[0]), at_ptp[1], at_ptp[2])
            if at_ptp[3]:
                s += '{:9.2f}{:9.2f}\n'.format(ALJ, ACRF)
            else:
                s += '{:9d}{:9d}\n'.format(0, 0)
        gf.write_block(bl, s)

    def __write_pert_interaction(self, gf, bl, pert_int_list, *args, **kwargs):
        s = str(len(pert_int_list)) + '\n'
        for pert_int in pert_int_list:
            for at in pert_int.atoms:
                s += at.write_gr_id
            for i in range(2):
                s += '{:>6}'.format(pert_int_list[pert_int][i])
            s += '\n'
        gf.write_block(bl, s)

    def __write_PERTATOMPAIR_v1(self, gf, bl, *args, **kwargs):
        s = str(len(self.ptp.pair_vdw)) + '\n'
        for at_p in self.ptp.pair_vdw:
            for at in at_p:
                s += at.write_gr_id
            for i in range(2):
                s += '{:>6}'.format(self.ptp.pair_vdw[at_p][i])
            s += '\n'
        gf.write_block(bl, s)

    def __write_PERTATOMPAIR_v2(self, gf, bl, *args, **kwargs):
        self.__write_pert_interaction(gf, bl, self.ptp.pair_vdw, *args, **kwargs)

    def __write_PERTBONDSTRETCH_v1(self, gf, bl, *args, **kwargs):
        self.__write_pert_interaction(gf, bl, self.ptp.bonds, *args, **kwargs)

    def __write_PERTBONDANGLE_v1(self, gf, bl, *args, **kwargs):
        self.__write_pert_interaction(gf, bl, self.ptp.angles, *args, **kwargs)

    def __write_PERTIMPROPERDIH_v1(self, gf, bl, *args, **kwargs):
        self.__write_pert_interaction(gf, bl, self.ptp.impropers, *args, **kwargs)

    def __write_PERTPROPERDIH_v1(self, gf, bl, *args, **kwargs):
        self.__write_pert_interaction(gf, bl, self.ptp.dihedrals, *args, **kwargs)

_ptpBlocksWriter_defs = {}
_ptpBlocksWriter_defs['_write_PERTATOMPARAM'] = '__write_PERTATOMPARAM_v1'
_ptpBlocksWriter_defs['_write_PERTATOMPAIR'] = '__write_PERTATOMPAIR_v1'
_ptpBlocksWriter_defs['_write_PERTBONDSTRETCH'] = '__write_PERTBONDSTRETCH_v1'
_ptpBlocksWriter_defs['_write_PERTBONDANGLE'] = '__write_PERTBONDANGLE_v1'
_ptpBlocksWriter_defs['_write_PERTIMPROPERDIH'] = '__write_PERTIMPROPERDIH_v1'
_ptpBlocksWriter_defs['_write_PERTPROPERDIH'] = '__write_PERTPROPERDIH_v1'

ptpBlocksWriter._add_defaults(_ptpBlocksWriter_defs, flag_set=1)
"""

# PTP / EDS reading and writing
class PTP_EDS_BlocksParser(GromosParser):
    def __MPERTATOM_v1(self, parse_from, bl_name):
        ff = self.get_intDB()
        N_EDS_atoms, N_states = (int(next(parse_from.block_split_fnc)), int(next(parse_from.block_split_fnc)))
        self.state_names = []
        for j in range(N_states):
            self.state_names.append(next(parse_from.block_split_fnc))
        for i in range(N_EDS_atoms):
            at = self.atoms[next(parse_from.block_split_fnc)]
            at_name = next(parse_from.block_split_fnc)
            assert at.name == at_name
            at.m_states = [at.m] * N_states
            at.p_ch_states = [None] * N_states
            at.a_type_states = [None] * N_states
            for j in range(N_states):
                at.a_type_states[j] = ff.a_type[next(parse_from.block_split_fnc)]
                at.p_ch_states[j] = float(next(parse_from.block_split_fnc))
            at.soft_core = float(next(parse_from.block_split_fnc)), float(next(parse_from.block_split_fnc))
        parse_from.end_block_split()

    def __write_MPERTATOM_v1(self, gf, bl, *args, **kwargs):
        DUM_type = self.get_DUM_type
        ALJ = kwargs.get('ALJ', 1.0)
        ACRF = kwargs.get('ACRF', 1.0)
        s = ''
        N_at_ptp = 0
        for at in self.get_atoms():
            if len(at.ND_a_type_states) != 1 or len(at.ND_pch_states) != 1:
                N_at_ptp += 1
                N_states = len(at.p_ch_states)
                s += at.write_gr_id + '{:>6}'.format(at.name)
                for i in range(N_states):
                    s += '{:>4}{:11.5f}'.format(get_id_string(at.a_type_states[i], flag_id='gr_id'), at.p_ch_states[i])
                s += '{:9.2f}{:9.2f}\n'.format(ALJ, ACRF)
        state_names = kwargs.get('state_names', [])
        if len(state_names)!=N_states:
            state_names = ['st_{:d}'.format(i) for i in range(1, N_states+1)]
        states_line = ' '*12 + '{:>15}' * N_states + '\n'
        s = '{:>5}{:>5}\n'.format(N_at_ptp, N_states) + states_line.format(*state_names) + s
        gf.write_block(bl, s)

    def write_EDS(self, f_path, **kwargs):
        temp_f = GromosFile(f_path, write=True)
        self.write_gromos_format(temp_f, 'TITLE')
        int_blocks2write = ('MPERTATOM', )
        self.write_gromos_format(temp_f, *int_blocks2write, **kwargs)


_PTP_EDS_BlocksParser_defs = {}
_PTP_EDS_BlocksParser_defs['_MPERTATOM_parser'] = '__MPERTATOM_v1'
_PTP_EDS_BlocksParser_defs['_write_MPERTATOM'] = '__write_MPERTATOM_v1'

PTP_EDS_BlocksParser._add_defaults(_PTP_EDS_BlocksParser_defs, flag_set=True)


## CNF reading and writing
class cnfBlocksParser(GromosParser):

    def __TIMESTEP_v1(self, parse_from, bl_name):
        self.time_step = (int(next(parse_from.block_split_fnc)), float(next(parse_from.block_split_fnc)))
        parse_from.end_block_split()

    @staticmethod
    def _get_pos_vel(line):
        temp = np.zeros(3)
        for i in range(3):
            temp[i] = float(line[24 + i * 15:39 + i * 15])
        return temp

    def __POSITION_v1(self, parse_from, bl_name):
        for i in parse_from.block_lines():
            if len(i) > 24:
                at = self.Atom()
                at.res_id = i[:5].strip()
                at.res_name = i[6:12].strip()
                at.name = i[12:18].strip()
                at.id = i[18:24].strip()
                at.coord = []
                at.coord = self._get_pos_vel(i)
                self.add2container(at, db_type=list, create = True)

    def __LATTICESHIFTS_v1(self, parse_from, bl_name):
        for at in self.atoms:
            at.latshift = []
            for i in range(3):
                at.latshift.append(next(parse_from.block_split_fnc))
        parse_from.end_block_split()

    def __VELOCITY_v1(self, parse_from, bl_name):
        c = 0
        for i in parse_from.block_lines():
            if len(i) > 24:
                at = self.atoms[c]
                at.vel = self._get_pos_vel(i)
                c += 1

    def __GENBOX_v1(self, parse_from, bl_name):
        self.box = self.Box()
        self.box.box_type = next(parse_from.block_split_fnc)
        for i in range(3):
            self.box.abc[i] = float(next(parse_from.block_split_fnc))
        for i in range(9):
            self.box.angles[i] = float(next(parse_from.block_split_fnc))
        parse_from.end_block_split()

_cnfBlocksParser_defs = {}
_cnfBlocksParser_defs['_TIMESTEP_parser'] = '__TIMESTEP_v1'
_cnfBlocksParser_defs['_POSITION_parser'] = '__POSITION_v1'
_cnfBlocksParser_defs['_LATTICESHIFTS_parser'] = '__LATTICESHIFTS_v1'
_cnfBlocksParser_defs['_VELOCITY_parser'] = '__VELOCITY_v1'
_cnfBlocksParser_defs['_GENBOX_parser'] = '__GENBOX_v1'
cnfBlocksParser._add_defaults(_cnfBlocksParser_defs, flag_set=True)


class cnfBlocksWriter(GromosWriter):

    @staticmethod
    def _get_pos_str(atoms):
        s = ''
        for at in atoms:
            s += '{:>5} {:<5} {:<6}{:>6}'.format(at.res_id, at.res_name, at.name, at.id)
            s += '{:15.9f}{:15.9f}{:15.9f}\n'.format(*at.coord)
        return s

    @staticmethod
    def _get_vel_str(atoms):
        s = ''
        for at in atoms:
            s += '{:>5} {:<5} {:<6}{:>6}'.format(at.res_id, at.res_name, at.name, at.id)
            if hasattr(at, 'vel'):
                s += '{:15.9f}{:15.9f}{:15.9f}\n'.format(*at.vel)
            else:
                s += '{:15.9f}{:15.9f}{:15.9f}\n'.format(0, 0, 0)
        return s

    def __write_TIMESTEP_v1(self, gf, bl, *args, **kwargs):
        s = '{:15d}{:15.9f}'.format(*self.time_step)
        gf.write_block(bl, s)

    def __write_POSITION_v1(self, gf, bl, *args, **kwargs):
        gf.write_block(bl, self._get_pos_str(self.atoms))

    def __write_VELOCITY_v1(self, gf, bl, *args, **kwargs):
        gf.write_block(bl, self._get_vel_str(self.atoms))

    def __write_POSRESSPEC_v1(self, gf, bl, *args, **kwargs):
        at2use = kwargs.get('at2use', [])
        if not at2use:
            at2discard = kwargs.get('at2discard', [])
            for at in self.atoms:
                if at.id not in at2discard:
                    at2use.append(at)
        else:
            if isinstance(at2use[0], int):
                at2use_cnf = []
                for at in at2use:
                    at2use_cnf.append(self.atoms[at])
                at2use = at2use_cnf
        s = self._get_pos_str(at2use)
        gf.write_block(bl, s)

    def __write_REFPOSITION_v1(self, gf, bl, *args, **kwargs):
        gf.write_block(bl, self._get_pos_str(self.atoms))

    def __write_GENBOX_v1(self, gf, bl, *args, **kwargs):
        s = '{:>5}\n'.format(self.box.box_type)
        temp = '{:15.9f}' * 3 + '\n'
        s += temp.format(*self.box.abc)
        for i in range(3):
            s += temp.format(*self.box.angles[i * 3:3 + i * 3])
        gf.write_block(bl, s)

_cnfBlocksWriter_defs = {}
_cnfBlocksWriter_defs['_write_POSITION'] = '__write_POSITION_v1'
_cnfBlocksWriter_defs['_write_VELOCITY'] = '__write_VELOCITY_v1'
_cnfBlocksWriter_defs['_write_TIMESTEP'] = '__write_TIMESTEP_v1'
_cnfBlocksWriter_defs['_write_POSRESSPEC'] = '__write_POSRESSPEC_v1'
_cnfBlocksWriter_defs['_write_REFPOSITION'] = '__write_REFPOSITION_v1'
_cnfBlocksWriter_defs['_write_GENBOX'] = '__write_GENBOX_v1'

cnfBlocksWriter._add_defaults(_cnfBlocksWriter_defs, flag_set=1)


#### CNF TRJ reading and writing

class TrjCnfBlocksParser(GromosParser):
    """parser for cnf trajectories"""
    def __just_for_readability(self):
        self.int_num_dtype = np.int32
        self.real_num_dtype = np.float32
        self.get_frame_dtype = None
        self._gr_add_frame = None

    @staticmethod
    def _get_pos_vel(line, real_num_dtype):
        temp = np.empty(3, dtype = real_num_dtype)
        for i, xyz in enumerate(split_gen(line)):
            temp[i] = xyz
        return temp

    def __TIMESTEP_v1(self, parse_from, bl_name, **kwargs):
        self.time_step = self.int_num_dtype(next(parse_from.block_split_fnc))
        self.time = self.time_real_num_dtype(next(parse_from.block_split_fnc))
        parse_from.end_block_split()

    def __POSITIONRED_v1(self, parse_from, bl_name, **kwargs):
        frame_coord = [line.split() for line in parse_from.block_lines()]
        self.coord = np.array(frame_coord, dtype=self.real_num_dtype)
        if self.N_atoms is None:
            self.N_atoms = len(frame_coord)
            self.get_frame_dtype()

    def __GENBOX_v1(self, parse_from, bl_name, **kwargs):
        self.box_type = np.int8(next(parse_from.block_split_fnc))
        temp_box =  np.empty((4,3), dtype=self.real_num_dtype)
        for i in range(4):
            for j in range(3):
                temp_box[i,j] = next(parse_from.block_split_fnc)
        parse_from.end_block_split()
        self.box = temp_box
        self._gr_add_frame()

_TrjCnfBlocksParser_defs = {}
_TrjCnfBlocksParser_defs['_TIMESTEP_parser'] = '__TIMESTEP_v1'
_TrjCnfBlocksParser_defs['_POSITIONRED_parser'] = '__POSITIONRED_v1'
_TrjCnfBlocksParser_defs['_GENBOX_parser'] = '__GENBOX_v1'
TrjCnfBlocksParser._add_defaults(_TrjCnfBlocksParser_defs, flag_set=True)


class GeneralCnfWriter(GromosWriter):
    def __just_for_readability(self):
        self._real_num_format = ''
        self._time_real_num_format = ''
        self.fr = None
        self._write_real_num_v = 1

    def __write_real_num(self, num, real_num_format):
        if num.shape:
            if len(num.shape) != 1:
                txt_format = (real_num_format * num.shape[1] + '\n') * num.shape[0]
                return txt_format.format(*num.flatten())
            else:
                txt_format = real_num_format * num.shape[0]
                return txt_format.format(*num)
        else:
            return real_num_format.format(num)

    def __write_real_num_v1(self, num):
        return self.__write_real_num(num, self._real_num_format)

    def __write_real_num_v2(self, num):
        return self.__write_real_num(num.astype(str), ' {:s}')

    def _get_atom_pos_str(self, atom_coord):
        fnc2write = getattr(self, '_GeneralCnfWriter__write_real_num_v' + str(self._write_real_num_v))
        return fnc2write(atom_coord) + '\n'

    def __write_TIMESTEP_v1(self, gf, bl, **kwargs):
        s = ('{:15d}' + self._time_real_num_format + '\n').format(self.fr[0][0], self.fr[1][0])
        gf.write_block(bl, s)

    def __write_GENBOX_v1(self, gf, bl, *args, **kwargs):
        s = '{:>5}\n'.format(self.fr['box_type'][0])
        fnc2write = getattr(self, '_GeneralCnfWriter__write_real_num_v' + str(self._write_real_num_v))
        s += fnc2write(self.fr['box'])
        gf.write_block(bl, s)

_GeneralCnfWriter_defs = {'_real_num_format':'{:15.9f}'}
_GeneralCnfWriter_defs = {'_time_real_num_format':'{:15.9f}'}
_GeneralCnfWriter_defs['_write_TIMESTEP'] = '__write_TIMESTEP_v1'
_GeneralCnfWriter_defs['_write_GENBOX'] = '__write_GENBOX_v1'
_GeneralCnfWriter_defs['_write_real_num_v'] = 1

GeneralCnfWriter._add_defaults(_GeneralCnfWriter_defs, flag_set=True)


class TrjCnfBlocksWriter(GeneralCnfWriter):
    def __just_for_readability(self):
        self.trj = None

    def write_trc_npz(self, f_path):
        if not f_path.endswith('.trc.npz'):
            if f_path.endswith('trc'):
                f_path += '.npz'
            else:
                f_path += '.trc.npz'
        np.savez_compressed(f_path, title = ''.join(self.sys_title[0].lines), trj = self.trj)

    def __write_POSITIONRED_v1(self, gf, bl, **kwargs):
        s = ''
        for at_coord in self.fr['coord']:
            s += self._get_atom_pos_str(at_coord)
        gf.write_block(bl, s)

    def write_trc(self, f_path = None, **kwargs):
        gs = self._get_grs_from_fpaht(f_path)
        self.write_gromos_format(gs, 'TITLE')
        int_blocks2write = ('TIMESTEP', 'POSITIONRED', 'GENBOX')
        for fr in self.trj:
            self.fr = fr
            self.write_gromos_format(gs, *int_blocks2write)
        if f_path and kwargs.get('flag_close', True):
            gs.f.close()

_TrjCnfBlocksWriter_defs = {}
_TrjCnfBlocksWriter_defs['_write_POSITIONRED'] = '__write_POSITIONRED_v1'

TrjCnfBlocksWriter._add_defaults(_TrjCnfBlocksWriter_defs, flag_set=True)
