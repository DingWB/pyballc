# This file was automatically generated by SWIG (http://www.swig.org).
# Version 4.0.2
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info < (2, 7, 0):
    raise RuntimeError("Python 2.7 or later required")

# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _pyballcools
else:
    import _pyballcools

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "thisown":
            self.this.own(value)
        elif name == "this":
            set(self, name, value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)


class SwigPyIterator(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _pyballcools.delete_SwigPyIterator

    def value(self):
        return _pyballcools.SwigPyIterator_value(self)

    def incr(self, n=1):
        return _pyballcools.SwigPyIterator_incr(self, n)

    def decr(self, n=1):
        return _pyballcools.SwigPyIterator_decr(self, n)

    def distance(self, x):
        return _pyballcools.SwigPyIterator_distance(self, x)

    def equal(self, x):
        return _pyballcools.SwigPyIterator_equal(self, x)

    def copy(self):
        return _pyballcools.SwigPyIterator_copy(self)

    def next(self):
        return _pyballcools.SwigPyIterator_next(self)

    def __next__(self):
        return _pyballcools.SwigPyIterator___next__(self)

    def previous(self):
        return _pyballcools.SwigPyIterator_previous(self)

    def advance(self, n):
        return _pyballcools.SwigPyIterator_advance(self, n)

    def __eq__(self, x):
        return _pyballcools.SwigPyIterator___eq__(self, x)

    def __ne__(self, x):
        return _pyballcools.SwigPyIterator___ne__(self, x)

    def __iadd__(self, n):
        return _pyballcools.SwigPyIterator___iadd__(self, n)

    def __isub__(self, n):
        return _pyballcools.SwigPyIterator___isub__(self, n)

    def __add__(self, n):
        return _pyballcools.SwigPyIterator___add__(self, n)

    def __sub__(self, *args):
        return _pyballcools.SwigPyIterator___sub__(self, *args)
    def __iter__(self):
        return self

# Register SwigPyIterator in _pyballcools:
_pyballcools.SwigPyIterator_swigregister(SwigPyIterator)

class AllC(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, file_path):
        _pyballcools.AllC_swiginit(self, _pyballcools.new_AllC(file_path))
    __swig_destroy__ = _pyballcools.delete_AllC

    def ReadLine(self):
        return _pyballcools.AllC_ReadLine(self)

# Register AllC in _pyballcools:
_pyballcools.AllC_swigregister(AllC)

class RefRecord(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    l_name = property(_pyballcools.RefRecord_l_name_get, _pyballcools.RefRecord_l_name_set)
    ref_name = property(_pyballcools.RefRecord_ref_name_get, _pyballcools.RefRecord_ref_name_set)
    ref_length = property(_pyballcools.RefRecord_ref_length_get, _pyballcools.RefRecord_ref_length_set)

    def __init__(self):
        _pyballcools.RefRecord_swiginit(self, _pyballcools.new_RefRecord())
    __swig_destroy__ = _pyballcools.delete_RefRecord

# Register RefRecord in _pyballcools:
_pyballcools.RefRecord_swigregister(RefRecord)

class BAllCHeader(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    magic = property(_pyballcools.BAllCHeader_magic_get, _pyballcools.BAllCHeader_magic_set)
    version = property(_pyballcools.BAllCHeader_version_get, _pyballcools.BAllCHeader_version_set)
    version_minor = property(_pyballcools.BAllCHeader_version_minor_get, _pyballcools.BAllCHeader_version_minor_set)
    sc = property(_pyballcools.BAllCHeader_sc_get, _pyballcools.BAllCHeader_sc_set)
    l_assembly = property(_pyballcools.BAllCHeader_l_assembly_get, _pyballcools.BAllCHeader_l_assembly_set)
    assembly_text = property(_pyballcools.BAllCHeader_assembly_text_get, _pyballcools.BAllCHeader_assembly_text_set)
    l_text = property(_pyballcools.BAllCHeader_l_text_get, _pyballcools.BAllCHeader_l_text_set)
    header_text = property(_pyballcools.BAllCHeader_header_text_get, _pyballcools.BAllCHeader_header_text_set)
    n_refs = property(_pyballcools.BAllCHeader_n_refs_get, _pyballcools.BAllCHeader_n_refs_set)
    refs = property(_pyballcools.BAllCHeader_refs_get, _pyballcools.BAllCHeader_refs_set)

    def __init__(self):
        _pyballcools.BAllCHeader_swiginit(self, _pyballcools.new_BAllCHeader())
    __swig_destroy__ = _pyballcools.delete_BAllCHeader

# Register BAllCHeader in _pyballcools:
_pyballcools.BAllCHeader_swigregister(BAllCHeader)
cvar = _pyballcools.cvar
BALLC_VERSION = cvar.BALLC_VERSION
BALLC_VERSION_MINOR = cvar.BALLC_VERSION_MINOR
BALLC_MAGIC = cvar.BALLC_MAGIC
BALLC_INDEX_MAGIC = cvar.BALLC_INDEX_MAGIC

class MCRecord(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    pos = property(_pyballcools.MCRecord_pos_get, _pyballcools.MCRecord_pos_set)
    ref_id = property(_pyballcools.MCRecord_ref_id_get, _pyballcools.MCRecord_ref_id_set)
    mc = property(_pyballcools.MCRecord_mc_get, _pyballcools.MCRecord_mc_set)
    cov = property(_pyballcools.MCRecord_cov_get, _pyballcools.MCRecord_cov_set)

    def __init__(self):
        _pyballcools.MCRecord_swiginit(self, _pyballcools.new_MCRecord())
    __swig_destroy__ = _pyballcools.delete_MCRecord

# Register MCRecord in _pyballcools:
_pyballcools.MCRecord_swigregister(MCRecord)

class MCRecord2(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    pos = property(_pyballcools.MCRecord2_pos_get, _pyballcools.MCRecord2_pos_set)
    chrom = property(_pyballcools.MCRecord2_chrom_get, _pyballcools.MCRecord2_chrom_set)
    mc = property(_pyballcools.MCRecord2_mc_get, _pyballcools.MCRecord2_mc_set)
    cov = property(_pyballcools.MCRecord2_cov_get, _pyballcools.MCRecord2_cov_set)

    def __init__(self):
        _pyballcools.MCRecord2_swiginit(self, _pyballcools.new_MCRecord2())
    __swig_destroy__ = _pyballcools.delete_MCRecord2

# Register MCRecord2 in _pyballcools:
_pyballcools.MCRecord2_swigregister(MCRecord2)

class BAllCFile(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    header = property(_pyballcools.BAllCFile_header_get, _pyballcools.BAllCFile_header_set)
    mc_records = property(_pyballcools.BAllCFile_mc_records_get, _pyballcools.BAllCFile_mc_records_set)

    def __init__(self):
        _pyballcools.BAllCFile_swiginit(self, _pyballcools.new_BAllCFile())
    __swig_destroy__ = _pyballcools.delete_BAllCFile

# Register BAllCFile in _pyballcools:
_pyballcools.BAllCFile_swigregister(BAllCFile)

class IndexHeader(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    magic = property(_pyballcools.IndexHeader_magic_get, _pyballcools.IndexHeader_magic_set)
    n_refs = property(_pyballcools.IndexHeader_n_refs_get, _pyballcools.IndexHeader_n_refs_set)
    refs = property(_pyballcools.IndexHeader_refs_get, _pyballcools.IndexHeader_refs_set)

    def __init__(self):
        _pyballcools.IndexHeader_swiginit(self, _pyballcools.new_IndexHeader())
    __swig_destroy__ = _pyballcools.delete_IndexHeader

# Register IndexHeader in _pyballcools:
_pyballcools.IndexHeader_swigregister(IndexHeader)

class IndexKey(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    ref_id = property(_pyballcools.IndexKey_ref_id_get, _pyballcools.IndexKey_ref_id_set)
    bin = property(_pyballcools.IndexKey_bin_get, _pyballcools.IndexKey_bin_set)

    def __lt__(self, other):
        return _pyballcools.IndexKey___lt__(self, other)

    def __eq__(self, other):
        return _pyballcools.IndexKey___eq__(self, other)

    def __init__(self):
        _pyballcools.IndexKey_swiginit(self, _pyballcools.new_IndexKey())
    __swig_destroy__ = _pyballcools.delete_IndexKey

# Register IndexKey in _pyballcools:
_pyballcools.IndexKey_swigregister(IndexKey)

class IndexEntry(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    key = property(_pyballcools.IndexEntry_key_get, _pyballcools.IndexEntry_key_set)
    chunk_start = property(_pyballcools.IndexEntry_chunk_start_get, _pyballcools.IndexEntry_chunk_start_set)
    chunk_end = property(_pyballcools.IndexEntry_chunk_end_get, _pyballcools.IndexEntry_chunk_end_set)

    def __lt__(self, other):
        return _pyballcools.IndexEntry___lt__(self, other)

    def __init__(self):
        _pyballcools.IndexEntry_swiginit(self, _pyballcools.new_IndexEntry())
    __swig_destroy__ = _pyballcools.delete_IndexEntry

# Register IndexEntry in _pyballcools:
_pyballcools.IndexEntry_swigregister(IndexEntry)

class BAllCIndexFile(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    header = property(_pyballcools.BAllCIndexFile_header_get, _pyballcools.BAllCIndexFile_header_set)
    n_entries = property(_pyballcools.BAllCIndexFile_n_entries_get, _pyballcools.BAllCIndexFile_n_entries_set)
    index = property(_pyballcools.BAllCIndexFile_index_get, _pyballcools.BAllCIndexFile_index_set)

    def __init__(self):
        _pyballcools.BAllCIndexFile_swiginit(self, _pyballcools.new_BAllCIndexFile())
    __swig_destroy__ = _pyballcools.delete_BAllCIndexFile

# Register BAllCIndexFile in _pyballcools:
_pyballcools.BAllCIndexFile_swigregister(BAllCIndexFile)

class BAllC(BAllCFile):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        _pyballcools.BAllC_swiginit(self, _pyballcools.new_BAllC(*args))
    __swig_destroy__ = _pyballcools.delete_BAllC

    def AllcLineToMcRecord(self, line):
        return _pyballcools.BAllC_AllcLineToMcRecord(self, line)

    def McRecordToLine(self, record):
        return _pyballcools.BAllC_McRecordToLine(self, record)

    def McRecordToMcRecord2(self, record):
        return _pyballcools.BAllC_McRecordToMcRecord2(self, record)

    def WriteMcRecord(self, record):
        return _pyballcools.BAllC_WriteMcRecord(self, record)

    def ReadMcRecord(self, record):
        return _pyballcools.BAllC_ReadMcRecord(self, record)

    def WriteHeader(self):
        return _pyballcools.BAllC_WriteHeader(self)

    def Seek(self, pos):
        return _pyballcools.BAllC_Seek(self, pos)

    def Tell(self):
        return _pyballcools.BAllC_Tell(self)

    def Close(self):
        return _pyballcools.BAllC_Close(self)
    ref_dict = property(_pyballcools.BAllC_ref_dict_get, _pyballcools.BAllC_ref_dict_set)

# Register BAllC in _pyballcools:
_pyballcools.BAllC_swigregister(BAllC)

class BAllCIndex(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, ballc_path):
        _pyballcools.BAllCIndex_swiginit(self, _pyballcools.new_BAllCIndex(ballc_path))
    __swig_destroy__ = _pyballcools.delete_BAllCIndex

    def BuildIndex(self):
        return _pyballcools.BAllCIndex_BuildIndex(self)

    def WriteIndex(self, override=False):
        return _pyballcools.BAllCIndex_WriteIndex(self, override)

    def QueryMcRecords_Iter(self, *args):
        return _pyballcools.BAllCIndex_QueryMcRecords_Iter(self, *args)
    GRANGE_END_MAX = _pyballcools.BAllCIndex_GRANGE_END_MAX

# Register BAllCIndex in _pyballcools:
_pyballcools.BAllCIndex_swigregister(BAllCIndex)

class MCRecordIterator(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, ballc, start_iter, end_iter, ref_id, start, end):
        _pyballcools.MCRecordIterator_swiginit(self, _pyballcools.new_MCRecordIterator(ballc, start_iter, end_iter, ref_id, start, end))

    def AdavanceToNextValid(self):
        return _pyballcools.MCRecordIterator_AdavanceToNextValid(self)

    def HasNext(self):
        return _pyballcools.MCRecordIterator_HasNext(self)

    def Next(self):
        return _pyballcools.MCRecordIterator_Next(self)

    def RefIdMatch(self, id1, id2):
        return _pyballcools.MCRecordIterator_RefIdMatch(self, id1, id2)
    ANY_REF_ID = _pyballcools.MCRecordIterator_ANY_REF_ID
    BAD_REF_ID = _pyballcools.MCRecordIterator_BAD_REF_ID
    __swig_destroy__ = _pyballcools.delete_MCRecordIterator

# Register MCRecordIterator in _pyballcools:
_pyballcools.MCRecordIterator_swigregister(MCRecordIterator)


def ExtractCMeta(fasta_path, cmeta_path):
    return _pyballcools.ExtractCMeta(fasta_path, cmeta_path)

def IndexCMeta(cmeta_path):
    return _pyballcools.IndexCMeta(cmeta_path)

def AllCToBallC(*args):
    return _pyballcools.AllCToBallC(*args)

def IndexBallc(ballc_path):
    return _pyballcools.IndexBallc(ballc_path)

def ExtractAllCMeta(fasta_path, cmeta_path):
    return _pyballcools.ExtractAllCMeta(fasta_path, cmeta_path)

def QueryBallc_Iter(ballc_path, range):
    return _pyballcools.QueryBallc_Iter(ballc_path, range)

def OutputMatched(*args):
    return _pyballcools.OutputMatched(*args)

def QueryBallcWithMeta_Iter(*args):
    return _pyballcools.QueryBallcWithMeta_Iter(*args)

def ViewBallc(*args):
    return _pyballcools.ViewBallc(*args)

def CheckBallc(ballc_path):
    return _pyballcools.CheckBallc(ballc_path)

def MergeBAllC(*args):
    return _pyballcools.MergeBAllC(*args)

def BAllCToAllC(*args):
    return _pyballcools.BAllCToAllC(*args)

def CombineCG(*args):
    return _pyballcools.CombineCG(*args)

def GetVersion():
    return _pyballcools.GetVersion()


