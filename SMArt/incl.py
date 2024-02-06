import sys
import os
import subprocess
import shutil
import copy
import math
import time
import gzip
from operator import itemgetter
import warnings
from warnings import warn
import traceback
import re
from re import finditer
from collections import OrderedDict, defaultdict, deque, Counter
from itertools import combinations, combinations_with_replacement, permutations, product
from bisect import bisect_left, bisect_right
from multiprocessing import Manager, Pool
import json
try:
    from frozendict import frozendict
except ImportError:
    frozendict = None
try:
    import pickle
except ImportError:
    pickle = None
try:
    import yaml
except ImportError:
    yaml = None
try:
    import numpy as np
except ImportError:
    np = None
try:
    import pandas as pd
except ImportError:
    pd = None
try:
    from matplotlib import pyplot as plt
except ImportError:
    plt = None
try:
    from scipy.optimize import minimize
except ImportError:
    minimize = None
try:
    from scipy import integrate
except ImportError:
    integrate = None

#import xml.etree.ElementTree as ET


class __VersionCompatibility():
    pass

if sys.version_info[0] < 3:
    class __VersionCompatibility(object):
        pass


#class __VersionCompatibility(object):pass

def test_time(f, *args, **kwargs):
    """
    calls f(*args, **kwargs) and times the execution (wall time)
    return:
        time
        res: return of the f(*args, **kwargs) call
    """
    t = time.time()
    res = f(*args, **kwargs)
    return time.time() - t, res


def split_gen(s, pattern=None):
    """the same as s.split() but it makes a generator out of it"""
    if not pattern:
        for m in finditer("\S+", s):
            yield m.group(0)
    else:
        for m in finditer(pattern, s):
            yield m.group(0)


class WARN(UserWarning):
    flag_traceback = 0

warnings.simplefilter('ignore', WARN)

def set_print_warnings(flag = 1):
    if flag:
        warnings.simplefilter('always', WARN)
    else:
        warnings.simplefilter('ignore', WARN)

def do_warn(m):
    if WARN.flag_traceback:
        temp_txt = ''
        temp_stack = traceback.format_stack()
        for i in temp_stack[:-1]:
            temp_txt+=i
        temp_warn = WARN('\n' + str(m) + '\n\ttraceback:\n' + temp_txt.strip())
        warn(temp_warn, stacklevel = 2)
    else:
        warn(m, WARN, stacklevel=2)

class Defaults(__VersionCompatibility):
    """
    takes care of default values:
        1) normal defaults (instance based), where each instance is independent
        2) class defaults where each instance inherit from the class, so change affects all instances
    """
    @classmethod
    def __transform_var(cls, var, flag_pvar = False):
        if isinstance(var, str):
            if var.startswith('__'):
                return '_' + cls.__name__ + var
            elif var.startswith('\__'):
                return var[1:]
            else:
                if not flag_pvar:
                    return var
                else:
                    return '_' + cls.__name__ + '__' + var
        else:
            return var

    @classmethod
    def __transform_defs(cls, defs, flag_pvar = False):
        new_dvalues = {}
        for var in defs:
            new_dvalues[cls.__transform_var(var, flag_pvar)] = cls.__transform_var(defs[var])
        return new_dvalues

    @classmethod
    def __add_class_defaults(cls, defs, default_dict, flag_set = False, flag_pvar = False):
        if not hasattr(cls, '_' + cls.__name__ + default_dict):
            setattr(cls, '_' + cls.__name__ + default_dict, dict())
        dvalues = getattr(cls, '_' + cls.__name__ + default_dict)
        new_dvalues = cls.__transform_defs(defs, flag_pvar)
        dvalues.update(new_dvalues)
        if flag_set:
            for var in new_dvalues:
                setattr(cls, var, new_dvalues[var])

    # normal defaults
    @classmethod
    def _add_defaults(cls, defs, flag_set = False):
        """add normal defaults, where defs is a dictionary - same format as __defaults: {var1:value, var2:value,...}"""
        cls.__add_class_defaults(defs, '__defaults', flag_set)

    @classmethod
    def __check_DefaultClass(cls):
        if hasattr(cls, "__mro__") and Defaults in cls.__mro__: return True
        return False

    @classmethod
    def __loop_base_Default_classes_mro(cls):
        for temp_base_class in cls.__mro__:
            if cls.__check_DefaultClass():
                yield temp_base_class

    @classmethod
    def __loop_base_Default_classes_defs_mro(cls):
        for kl in cls.__loop_base_Default_classes_mro():
            yield getattr(kl, '_' + kl.__name__ + '__defaults', {})

    @classmethod
    def get_defaults(cls, *vars):
        """get normal (instance-based) defaults - list all options"""
        dvalues = {}
        if not vars:
            for dd in cls.__loop_base_Default_classes_defs_mro():
                for var in dd:
                    if var not in dvalues:
                        dvalues[var] = []
                    dvalues[var].append(dd[var])
        else:
            for dd in cls.__loop_base_Default_classes_defs_mro():
                for var in dd:
                    if var in vars:
                        if var not in dvalues:dvalues[var] = []
                        dvalues[var].append(dd[var])
        for var in vars:
            if var not in dvalues:
                do_warn('variable {:} not in defaults'.format(var))
        return dvalues

    def set_defaults(self, *vars, **kwargs):
        """set normal (instance-based) defaults"""
        version = kwargs.get('version', 0)
        dd = self.get_defaults(*vars)
        for var in dd:
            setattr(self, var, dd[var][version])

    def get_current_dvalues(self):
        """get current values of default variables (normal, instance-based)"""
        current_dvalues = {}
        dvalues = self.get_defaults()
        for var in dvalues:
            if hasattr(self, var):
                current_dvalues[var] = getattr(self, var)
        return current_dvalues

    def reset_dvalues(self, temp_dvalues):
        """reset values of default variables from some state obtained using get_current_dvalues function"""
        dvalues = self.get_defaults()
        for var in dvalues:
            if var in temp_dvalues:
                setattr(self, var, temp_dvalues[var])
            else:
                if hasattr(self, var):
                    try:
                        delattr(self, var)
                    except:
                        do_warn('not able to remove this variable: '+var)

    # class defaults
    @classmethod
    def _add_class_defaults(cls, defs, flag_set = False):
        """add class defaults, where defs is a dictionary
        same format as __class_defaults: {__var1:value, __var2:value,...}"""
        cls.__add_class_defaults(defs, '__class_defaults', flag_set, flag_pvar = True)

    @classmethod
    def _get_class_defaults(cls, *vars, **kwargs):
        """
        get class defaults - list all options
        :param vars:
        :param kwargs:
            flag_do_warn = True - prints a warning if the class doesn't have defaults
        :return:
        """
        d_values = {}
        d_var = '_' + cls.__name__ + '__class_defaults'
        if hasattr(cls, '_' + cls.__name__ + '__class_defaults'):
            d_dict = getattr(cls,d_var)
            if not vars:
                vars = d_dict.keys()
            for var in vars:
                if not var.startswith('_' + cls.__name__ + '__'):
                    var = '_' + cls.__name__ + '__' + var
                if var in d_dict:
                    d_values[var] = d_dict[var]
                else:
                    do_warn('variable {:} not in defaults'.format(var))
        else:
            if kwargs.get('flag_do_warn', True):
                do_warn("class {:} doesn't have defaults".format(cls.__name__))
        return d_values

    @classmethod
    def _set_class_defaults(cls, *vars, **kwargs):
        """
        set class defaults
        :param vars:
        :param kwargs:
            flag_do_warn = True - prints a warning if the class doesn't have defaults
        :return:
        """
        d_values = cls._get_class_defaults(*vars, **kwargs)
        for var in d_values:
            setattr(cls, var, d_values[var])

    @classmethod
    def _get_current_class_dvalues(cls, **kwargs):
        """
        get current values of default class variables
        :param kwargs:
            flag_do_warn = True - prints a warning if the class doesn't have defaults
        :return:
        """
        curr_values = {}
        d_values = cls._get_class_defaults(**kwargs)
        for var in d_values:
            if hasattr(cls, var):
                curr_values[var] = getattr(cls, var)
        return curr_values

    @classmethod
    def _reset_class_dvalues(cls, temp_dvalues):
        """reset values of default class variables from some state obtained using _get_current_class_dvalues function"""
        d_values = cls._get_class_defaults()
        for var in d_values:
            if var in temp_dvalues:
                setattr(cls, var, temp_dvalues[var])
            else:
                if hasattr(cls, var):
                    try:
                        delattr(cls, var)
                    except:
                        do_warn('not able to remove this variable: '+var)

class FileHandler(Defaults):
    def __init__(self, f_path, write=False, bck_num=False, pre=False, suf=False, flag_bck=True, **kwargs):
        """creates or reads a file at f_path
        if bck_num given, reads the back-up file
        if write given, will back up a file if the file already exists
        pre and suf are used for naming back-up files: bck_file_name = pre + file_name + suf + str(bck_num)"""
        self.set_defaults()
        if pre:
            self.__change_pre(pre)
        if suf:
            self.__change_suf(suf)
        self.__new_f(f_path, write, bck_num, flag_bck, **kwargs)
        self.f_path = os.path.abspath(f_path)

    __defaults = {'pre':'#bck_', 'suf':'_'}

    def __change_pre(self, pre):
        """change prefix for backup files"""
        self.pre = pre

    # change suf
    def __change_suf(self, suf):
        """change sufix for backup files"""
        self.suf = suf

    def get_f(self):
        """returns the file"""
        return self.f

    def __new_f(self, f_path, write=False, bck_num=False, flag_bck=True, fnc2open = open, **kwargs):
        assert isinstance(f_path, str)
        if write and bck_num:
            raise Exception("write and bck_num mutually exclusive")
        if write:
            if os.path.isfile(f_path) and flag_bck:
                self.__bck_file(f_path)
            try:
                mode = kwargs.get('mode', 'w')
                if not mode.startswith('w'):
                    raise Exception('open mode has to start with "w" - otherwise set write = False')
                self.f = fnc2open(f_path, mode)
            except IOError:
                raise Exception("currently in: " + os.getcwd() + "\ncannot open: " + f_path)
        else:
            if bck_num:
                f_path = self.__bck_file_name(os.path.split(f_path)[0], os.path.split(f_path)[1], bck_num)
            try:
                mode = kwargs.get('mode', 'r')
                if not mode.startswith('r'):
                    raise Exception('open mode has to start with "r" - otherwise set write = True')
                self.f = fnc2open(f_path, mode)
            except IOError:
                raise Exception("currently in: " + os.getcwd() + "\ncannot open: " + f_path)

    def __bck_file(self, f_path):
        """back up a file - called automatically when a file gets open with write=True"""
        fd, f_name = os.path.split(f_path)
        c = 1
        temp_f = self.__bck_file_name(fd, f_name, c)
        while os.path.isfile(temp_f):
            c += 1
            temp_f = self.__bck_file_name(fd, f_name, c)
        shutil.move(f_path, temp_f)

    def __bck_file_name(self, fd, f_name, c):
        """returns a name for a backup file - called from bck_file"""
        return os.path.join(fd, self.pre + f_name + self.suf + str(c))


def _copy_file(fname, fd, cp_fname = None):
    if cp_fname is None:
        cp_fname = os.path.split(fname)[1]
    shutil.copy(fname, fd + '/' + cp_fname)
    return cp_fname


def _parse_line(self, l):
    """removes comments"""
    if isinstance(l, bytes):
        l = l.decode("utf-8")
    for comm in self._comm_indicators:
        if comm in l:
            l = l[:l.find(comm)]
    return l


class FileStream(FileHandler):
    parse_line = _parse_line

    def lines(self):
        """generator yielding line by line (comments removed)"""
        for l in self.f:
            temp_l = self.parse_line(l)
            if temp_l:
                yield (temp_l)

    def reset_stream(self):
        self.f = open(self.f_path)

    def _remove_f(self):
        if getattr(self, 'f'):
            self.f = None


class StringStream(Defaults):
    f = sys.stdout
    parse_line = _parse_line

    def __init__(self, s = None):
        if s is None:
            s = ''
        if isinstance(s, str):
            self.s = s
            self.__s = str(s)
            self.f_path = None
        else:
            raise TypeError('string expected - provided ' + type(s).__name__)

    def reset_stream(self):
        self.s = str(self.__s)

    def lines(self):
        """generator yielding line by line (comments removed) - from string"""
        while self.s:
            temp = self.s.split('\n', 1)
            temp_l = temp.pop(0)
            if temp:
                self.s = temp.pop()
            else:
                self.s = ''
            yield(self.parse_line(temp_l))

    def _remove_f(self):
        return

class ContainerItem:
    def __just_for_readability(self):
        self.container2write = None # this defines how to call a container
        self.id = None # id of the item
        self.name # name of the item
        self._name_var = None

    @property
    def __prf(self):
        txt = str(self.id)
        if hasattr(self, 'name'):
            name = str(self.name)
        else:
            _name_var = getattr(self, '_name_var')
            if _name_var:
                name = getattr(self, _name_var)
        if name:
            txt += ' ' + name
        return txt

    def __str__(self):
        return self.__prf

    def __repr__(self):
        return self.__prf


class GeneralContainer:
    ContainerItem = ContainerItem

    def get_container(self, container_name, create=False, db_type=OrderedDict, **kwargs):
        """gets the container by the name
        kwargs: 
            flag_item - searches for the container name from the item (e.g. flag_item = Bond())
            flag_class - searches for the container name from the class (e.g. flag_item = Bond)
                flag_item and flag_class finds the container name that is usually defined from container2write
            pre_containers - list of pre_containers (e.g. top.ff with cont_name = bonds gives top.ff.bonds)
            cont_pref - prefix for container name
            cont_suf - sufix for container name
            allow_not_found - if container is not found, Exception is raised, 
                if this flag is set to True, this prevents the Exception
        """
        if kwargs.get('flag_item', False):
            container_name = self.__get_container_name(container_name, **kwargs)
        if kwargs.get('flag_class', False):
            container_name = self.__get_container_name_class(container_name, **kwargs)
        if type(container_name) is not list:
            if type(container_name) is str:
                container_name = [container_name]
            else:
                container_name = list(container_name)
        container_names = kwargs.get('pre_containers', [])
        container_names.extend(container_name)
        return self.__get_container_recursive(container_names, create, db_type, **kwargs)

    def __update_containers(self, container_name):
        if not hasattr(self, '_containers'):
            self._containers = []
        if container_name not in self._containers:
            self._containers.append(container_name)

    def __get_container_single(self, container_name, create=False, db_type=OrderedDict, **kwargs):
        if not hasattr(self, container_name):
            if create:
                setattr(self, container_name, db_type())
                self.__update_containers(container_name)
            elif kwargs.get('allow_not_found'):
                return
            else:
                raise Exception('container ' + container_name + ' not defined; set create=True')
        return getattr(self, container_name)

    def __get_container_group(self, container_name, create=False, **kwargs):
        if not hasattr(self, container_name):
            if create:
                setattr(self, container_name, GeneralContainer())
                self.__update_containers(container_name)
            elif kwargs.get('allow_not_found'):
                return
            else:
                raise Exception('container group ' + container_name + ' not defined; set create=True')
        return getattr(self, container_name)

    def __get_container_recursive(self, container_names, create=False, db_type=OrderedDict, **kwargs):
        container_name = container_names.pop(0)
        if container_names:
            db = self.__get_container_group(container_name, create, **kwargs)
            return db.__get_container_recursive(container_names, create, db_type, **kwargs)
        else:
            return self.__get_container_single(container_name, create, db_type, **kwargs)

    @staticmethod
    def __get_pref_suf_container(item_class, cont_pref = None, cont_suf = None, **kwargs):
        if cont_pref is None:
            if not hasattr(item_class, 'cont_pref'):
                cont_pref = ''
            else:
                cont_pref = getattr(item_class, 'cont_pref')
        if cont_suf is None:
            if not hasattr(item_class, 'cont_suf'):
                cont_suf = ''
            else:
                cont_suf = getattr(item_class, 'cont_suf')
        return cont_pref, cont_suf

    def __get_container_name_from_container2write(self, item_class, **kwargs):
        cont_pref, cont_suf = self.__get_pref_suf_container(item_class, **kwargs)
        return cont_pref + getattr(item_class, 'container2write') + cont_suf

    def __get_container_name(self, item, **kwargs):
        if hasattr(item, 'container2write'):
            return self.__get_container_name_from_container2write(item, **kwargs)
        return 'db_' + item.__class__.__name__

    def __get_container_name_class(self, klass, **kwargs):
        if hasattr(klass, 'container2write'):
            return self.__get_container_name_from_container2write(klass, **kwargs)
        return 'db_' + klass.__name__

    @staticmethod
    def __get_pref_suf_item_id(klass):
        if not hasattr(klass, 'id_pref'):
            id_pref = ''
        else:
            id_pref = getattr(klass, 'id_pref')
        if not hasattr(klass, 'id_suf'):
            id_suf = ''
        else:
            id_suf = getattr(klass, 'id_suf')
        return id_pref, id_suf

    def add2container(self, item, replace=None, create=False, container_name=None, db_type=OrderedDict, **kwargs):
        """adds an item to a container, container name comes from the item, but can be passed as a variable
        kwargs: list_index - if the container is a list, the item will get inserted to the given index
                item_id - can be passed if the item doesnt have its own it
                replace = (-1, 0, 1) -1 used to ignore if an item is already in the container, 0 raises an error
                flag_find_next_id"""
        if not container_name:
            container_name = self.__get_container_name(item)
        db = self.get_container(container_name, create, db_type, **kwargs)
        from_segment = kwargs.get('from_segment', None)
        if from_segment:
            item.from_segment = from_segment
        if type(db) is list:
            if item in db:
                if replace == 1:
                    db.remove(item)
                elif replace==-1:
                    return
                else:
                    raise Exception(
                        'item already in the list' + str(item.id) + '; add replace = True or replace = -1 (ignore)')
            ind = kwargs.get('list_index', len(db))
            db.insert(ind, item)
        elif set in db.__class__.__mro__:
            db.add(item)
        else:
            if not hasattr(item, 'id') and 'item_id' not in kwargs:
                if kwargs.get('flag_find_next_id', True):
                    item.id = self._find_next_code(db, **kwargs)
                else:
                    raise Exception('give item id to add it to a container (either as item.id or item_id variable of the'
                                    ' function) or use flag_find_next_id')
            if hasattr(item, 'id'):
                item_id = item.id
            else:
                item_id = kwargs['item_id']
            if item_id is None:
                if kwargs.get('flag_find_next_id', True):
                    item.id = self._find_next_code(db, **kwargs)
                    item_id = item.id
                else:
                    raise Exception('item id is None.... set it to something else or use flag_find_next_id')
            if item_id in db:
                if replace==1:
                    db[item_id] = item
                elif replace==-1:
                    return
                else:
                    raise Exception(str(item_id) + ' already defined in the database ' + str(
                                    container_name) + '; add replace = True or replace = -1 (ignore)')
            else:
                db[item_id] = item

    def get_item(self, item_id, klass, create=False, create_container=False, **kwargs):
        db = self.get_container(klass, flag_class=True, create=create_container, **kwargs)
        if dict in db.__class__.__mro__:
            if item_id is not None and item_id in db:
                return db[item_id]
        if create:
            if dict in db.__class__.__mro__:
                if item_id is None:
                    id_pref, id_suf = self.__get_pref_suf_item_id(klass)
                    item_id = self._find_next_code(db, pref=id_pref, suf=id_suf, **kwargs)
                item = klass(**kwargs)
                item.id = item_id
            else:
                item = klass(**kwargs)
            self.add2container(item, **kwargs)
            return item
        if kwargs.get('allow_not_found'):
            return
        raise Exception('item_id not defined: ' + str(item_id))

    @staticmethod
    def get_code(c, pref = '', suf = '', **kwargs):
        return pref + str(c) + suf

    @staticmethod
    def get_num_from_code(code, pref = '', suf = '', **kwargs):
        pos2 = -len(suf)
        if not pos2:
            pos2=None
        return int(code[len(pref):pos2])

    def _find_next_code(self, container, c=1, pref=None, suf=None, **kwargs):
        id_pref, id_suf = '', ''
        klass = kwargs.get('klass')
        if klass:
            id_pref, id_suf = self.__get_pref_suf_item_id(klass)
        if pref is None:pref = id_pref
        if suf is None:suf = id_suf
        temp_code = self.get_code(c, pref, suf)
        flag_type = kwargs.get('id_type')
        if flag_type:
            temp_code = flag_type(temp_code)
        while temp_code in container:
            c += 1
            temp_code = self.get_code(c, pref, suf)
            if flag_type:
                temp_code = flag_type(temp_code)
        return temp_code

    def find_next_code(self, container_name, c=1, pref=None, suf=None, create=False, **kwargs):
        """
        finds the next available code in a dict (of e.g. bonds)
        :param kwargs:
            id_type - defines the type
        :return:
        """
        db = self.get_container(container_name, create, **kwargs)
        return self._find_next_code(db, c, pref, suf, **kwargs)

    @staticmethod
    def __renumber_list(temp_cont, c_init, attrib, flag_id_map, **kwargs):
        id_map = {}
        c = c_init
        for temp_item in temp_cont:
            new_num = kwargs['id_type'](c)
            if flag_id_map:
                id_map[getattr(temp_item, attrib)] = new_num
            setattr(temp_item, attrib, new_num)
            c+=1
        return id_map

    @staticmethod
    def __renumber_ordereddict(temp_cont, c_init, attrib, flag_id_map, **kwargs):
        id_map = {}
        c = c_init
        temp_list_items = []
        for temp_item_id in list(temp_cont):
            temp_item = temp_cont.pop(temp_item_id)
            new_num = kwargs['id_type'](c)
            if flag_id_map:
                id_map[getattr(temp_item, attrib)] = new_num
            setattr(temp_item, attrib, new_num)
            temp_list_items.append(temp_item)
            c+=1
        for temp_item in temp_list_items:
            temp_cont[temp_item.id] = temp_item
        return id_map

    def renumber_container(self, container_name, c_init = 1, attrib = 'id', flag_id_map = False, **kwargs):
        """
        :param container_name:
        :param c_init:
        :param attrib:
        :param flag_id_map:
        :param kwargs:
            id_type - str or int
        :return:
        """
        db = self.get_container(container_name, **kwargs)
        kwargs['id_type'] = kwargs.get('id_type', str)
        if OrderedDict in db.__class__.__mro__:
            return self.__renumber_ordereddict(db, c_init, attrib, flag_id_map, **kwargs)
        elif list in db.__class__.__mro__ or tuple in db.__class__.__mro__:
            return self.__renumber_list(db, c_init, attrib, flag_id_map, **kwargs)
        else:
            raise Exception('type(container) not in (list, tuple, OrderedDict)...')

    def __sort_dict(self, db, attrib, pref, suf, **kwargs):
        if attrib is None:
            temp_list_ids = [(self.get_num_from_code(code, pref, suf), code) for code in db]
        else:
            temp_list_ids = [(self.get_num_from_code(getattr(db[item_id], attrib), pref, suf), item_id) for item_id in db]
        temp_list_ids.sort(reverse = kwargs.get('reverse', False))
        for c, code in temp_list_ids:
            temp_item = db.pop(code)
            db[code] = temp_item

    def __sort_list(self, db, attrib, pref, suf, **kwargs):
        db.sort(key=lambda item: self.get_num_from_code(getattr(item, attrib), pref, suf), reverse = kwargs.get('reverse', False))

    def sort_container(self, container_name, attrib = 'id', pref='', suf='', **kwargs):
        db = self.get_container(container_name, **kwargs)
        if OrderedDict in db.__class__.__mro__:
            self.__sort_dict(db, attrib, pref, suf)
        elif list in db.__class__.__mro__:
            self.__sort_list(db, attrib, pref, suf)
        else:
            raise Exception('type(container) not in (list, OrderedDict)...')

def get_id_string(item, flag_id = 'id', default_id = 'id'):
    if hasattr(item, flag_id):
        return getattr(item, flag_id)
    if flag_id!=default_id and hasattr(item, default_id):
        return getattr(item, default_id)
    return item

mode_map = {'pickle':'b'}


class DataDumping:
    """class providing pickle of yaml dump functionality"""
    def __pickle_dump(self, temp_f):
        pickle.dump(self, temp_f.f)

    def __yaml_dump(self, temp_f):
        yaml.dump(self, temp_f.f)

    def __remove_files(self, **kwargs):
        extra_parsing_sources = kwargs.get('extra_parsing_sources', ['ff'])
        if hasattr(self, 'parsed_from'):
            for parsing_source in self.parsed_from:
                if FileHandler in parsing_source.__class__.__mro__ and parsing_source.f is not None:
                    parsing_source.f = None
        for temp_attr in extra_parsing_sources:
            temp_parsed_source = getattr(self, temp_attr, None)
            if temp_parsed_source and hasattr(temp_parsed_source, '_DataDumping__remove_files'):
                temp_parsed_source._DataDumping__remove_files(**kwargs)

    def dump(self, f_path, dump_format = 'pickle', **kwargs):
        """
        dumps data in a file
        :param f_path:
        :param dump_format: ('pickle', 'yaml')
        :param kwargs:
            mode: same as file open mode
            prep4dump - additional function to call before dumping
            post_dump_load - additional function to call after dumping (usually opposite from prep4dump)
            other kwargs from SMArt.incl.FileHandler
        :return:
        """
        if dump_format in mode_map and 'mode' not in kwargs:
            kwargs['mode'] = 'w' + mode_map[dump_format]
        fnc_map = {'pickle':self.__pickle_dump,'yaml':self.__yaml_dump}
        if dump_format not in fnc_map:
            raise Exception("dump_format not in", fnc_map.keys())
        if kwargs.get('prep4dump'):
            prep4dump_fnc = getattr(self, kwargs.get('prep4dump'), None)
            if prep4dump_fnc:
                prep4dump_fnc(**kwargs)
        #self.__remove_files()
        temp_f = FileHandler(f_path, write = True, **kwargs)
        fnc_map[dump_format](temp_f)
        temp_f.f.close()
        if kwargs.get('post_dump_load'):
            post_dump_load_fnc = getattr(self, kwargs.get('post_dump_load'), None)
            if post_dump_load_fnc:
                post_dump_load_fnc(**kwargs)

    @staticmethod
    def __pickle_load(temp_f):
        return pickle.load(temp_f.f)

    @staticmethod
    def __yaml_load(temp_f):
        return yaml.load(temp_f.f, Loader=yaml.Loader)

    @classmethod
    def load(cls, f_path, load_format = 'pickle', **kwargs):
        """
        loads data from a file
        :param f_path:
        :param load_format: ('pickle', 'yaml')
        :param kwargs:
            mode: same as file open mode
            post_dump_load - additional function to call after loading
            other kwargs from SMArt.incl.FileHandler
        :return: data
        """
        if load_format in mode_map and 'mode' not in kwargs:
            kwargs['mode'] = 'r' + mode_map[load_format]
        fnc_map = {'pickle': cls.__pickle_load, 'yaml': cls.__yaml_load}
        if load_format not in fnc_map:
            raise Exception("dump_format not in", fnc_map.keys())
        temp_f = FileHandler(f_path, write=False, **kwargs)
        loaded_data = fnc_map[load_format](temp_f)
        temp_f.f.close()
        if kwargs.get('post_dump_load'):
            post_dump_load_fnc = getattr(loaded_data, kwargs.get('post_dump_load'), None)
            if post_dump_load_fnc:
                post_dump_load_fnc(**kwargs)
        return loaded_data

import argparse
class ArgParser(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_line):
        return arg_line.split()