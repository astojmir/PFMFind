# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.
import _FS
def _swig_setattr(self,class_type,name,value):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    self.__dict__[name] = value

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class seqn(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, seqn, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, seqn, name)
    def __init__(self,*args):
        _swig_setattr(self, seqn, 'this', apply(_FS.new_seqn,args))
        _swig_setattr(self, seqn, 'thisown', 1)
    def __del__(self, destroy= _FS.delete_seqn):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __str__(*args): return apply(_FS.seqn___str__,args)
    __swig_getmethods__["len"] = _FS.seqn_len_get
    if _newclass:len = property(_FS.seqn_len_get)
    __swig_getmethods__["seq"] = _FS.seqn_seq_get
    if _newclass:seq = property(_FS.seqn_seq_get)
    __swig_getmethods__["defline"] = _FS.seqn_defline_get
    if _newclass:defline = property(_FS.seqn_defline_get)
    def __repr__(self):
        return "<C seqn instance at %s>" % (self.this,)

class seqnPtr(seqn):
    def __init__(self,this):
        _swig_setattr(self, seqn, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, seqn, 'thisown', 0)
        _swig_setattr(self, seqn,self.__class__,seqn)
_FS.seqn_swigregister(seqnPtr)

class db(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, db, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, db, name)
    __swig_getmethods__["db_name"] = _FS.db_db_name_get
    if _newclass:db_name = property(_FS.db_db_name_get)
    __swig_getmethods__["no_seq"] = _FS.db_no_seq_get
    if _newclass:no_seq = property(_FS.db_no_seq_get)
    def __init__(self,*args):
        _swig_setattr(self, db, 'this', apply(_FS.new_db,args))
        _swig_setattr(self, db, 'thisown', 1)
    def __del__(self, destroy= _FS.delete_db):
        try:
            if self.thisown: destroy(self)
        except: pass
    def __str__(*args): return apply(_FS.db___str__,args)
    def get_seq(*args): return apply(_FS.db_get_seq,args)
    def init_frags(*args): return apply(_FS.db_init_frags,args)
    def clear_frags(*args): return apply(_FS.db_clear_frags,args)
    def get_nofrags(*args): return apply(_FS.db_get_nofrags,args)
    def get_frag(*args): return apply(_FS.db_get_frag,args)
    def init_Ffrags(*args): return apply(_FS.db_init_Ffrags,args)
    def count_Ffrags(*args): return apply(_FS.db_count_Ffrags,args)
    def get_Ffrag(*args): return apply(_FS.db_get_Ffrag,args)
    def __repr__(self):
        return "<C db instance at %s>" % (self.this,)

class dbPtr(db):
    def __init__(self,this):
        _swig_setattr(self, db, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, db, 'thisown', 0)
        _swig_setattr(self, db,self.__class__,db)
_FS.db_swigregister(dbPtr)

class ptable(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ptable, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ptable, name)
    def __init__(self,*args):
        _swig_setattr(self, ptable, 'this', apply(_FS.new_ptable,args))
        _swig_setattr(self, ptable, 'thisown', 1)
    def __del__(self, destroy= _FS.delete_ptable):
        try:
            if self.thisown: destroy(self)
        except: pass
    def print_table(*args): return apply(_FS.ptable_print_table,args)
    def __str__(*args): return apply(_FS.ptable___str__,args)
    __swig_getmethods__["no_pttn"] = _FS.ptable_no_pttn_get
    if _newclass:no_pttn = property(_FS.ptable_no_pttn_get)
    def get_pttn(*args): return apply(_FS.ptable_get_pttn,args)
    def get_posn(*args): return apply(_FS.ptable_get_posn,args)
    def get_pttn_size(*args): return apply(_FS.ptable_get_pttn_size,args)
    def get_poffset(*args): return apply(_FS.ptable_get_poffset,args)
    def get_letter(*args): return apply(_FS.ptable_get_letter,args)
    def check_seqn(*args): return apply(_FS.ptable_check_seqn,args)
    def write(*args): return apply(_FS.ptable_write,args)
    def seqn2reduced(*args): return apply(_FS.ptable_seqn2reduced,args)
    def __repr__(self):
        return "<C ptable instance at %s>" % (self.this,)

class ptablePtr(ptable):
    def __init__(self,this):
        _swig_setattr(self, ptable, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ptable, 'thisown', 0)
        _swig_setattr(self, ptable,self.__class__,ptable)
_FS.ptable_swigregister(ptablePtr)

ptable_read = _FS.ptable_read

class fgen(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, fgen, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, fgen, name)
    __swig_getmethods__["max_len"] = _FS.fgen_max_len_get
    if _newclass:max_len = property(_FS.fgen_max_len_get)
    def __init__(self,*args):
        _swig_setattr(self, fgen, 'this', apply(_FS.new_fgen,args))
        _swig_setattr(self, fgen, 'thisown', 1)
    def __del__(self, destroy= _FS.delete_fgen):
        try:
            if self.thisown: destroy(self)
        except: pass
    def rand_seq(*args): return apply(_FS.fgen_rand_seq,args)
    def __repr__(self):
        return "<C fgen instance at %s>" % (self.this,)

class fgenPtr(fgen):
    def __init__(self,this):
        _swig_setattr(self, fgen, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, fgen, 'thisown', 0)
        _swig_setattr(self, fgen,self.__class__,fgen)
_FS.fgen_swigregister(fgenPtr)

class smatrix(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, smatrix, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, smatrix, name)
    __swig_getmethods__["filename"] = _FS.smatrix_filename_get
    if _newclass:filename = property(_FS.smatrix_filename_get)
    __swig_getmethods__["similarity_flag"] = _FS.smatrix_similarity_flag_get
    if _newclass:similarity_flag = property(_FS.smatrix_similarity_flag_get)
    def __init__(self,*args):
        _swig_setattr(self, smatrix, 'this', apply(_FS.new_smatrix,args))
        _swig_setattr(self, smatrix, 'thisown', 1)
    def __del__(self, destroy= _FS.delete_smatrix):
        try:
            if self.thisown: destroy(self)
        except: pass
    def set_M(*args): return apply(_FS.smatrix_set_M,args)
    def set_SS(*args): return apply(_FS.smatrix_set_SS,args)
    def get_M(*args): return apply(_FS.smatrix_get_M,args)
    def get_pM(*args): return apply(_FS.smatrix_get_pM,args)
    def get_pMc(*args): return apply(_FS.smatrix_get_pMc,args)
    def get_SS(*args): return apply(_FS.smatrix_get_SS,args)
    def score(*args): return apply(_FS.smatrix_score,args)
    def print_matrix(*args): return apply(_FS.smatrix_print_matrix,args)
    def S2Dmax(*args): return apply(_FS.smatrix_S2Dmax,args)
    def S2Davg(*args): return apply(_FS.smatrix_S2Davg,args)
    def S2Dquasi(*args): return apply(_FS.smatrix_S2Dquasi,args)
    def __repr__(self):
        return "<C smatrix instance at %s>" % (self.this,)

class smatrixPtr(smatrix):
    def __init__(self,this):
        _swig_setattr(self, smatrix, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, smatrix, 'thisown', 0)
        _swig_setattr(self, smatrix,self.__class__,smatrix)
_FS.smatrix_swigregister(smatrixPtr)

class hit(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, hit, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, hit, name)
    __swig_getmethods__["sequence_id"] = _FS.hit_sequence_id_get
    if _newclass:sequence_id = property(_FS.hit_sequence_id_get)
    __swig_getmethods__["sequence_from"] = _FS.hit_sequence_from_get
    if _newclass:sequence_from = property(_FS.hit_sequence_from_get)
    __swig_getmethods__["rejected"] = _FS.hit_rejected_get
    if _newclass:rejected = property(_FS.hit_rejected_get)
    __swig_getmethods__["value"] = _FS.hit_value_get
    if _newclass:value = property(_FS.hit_value_get)
    __swig_getmethods__["pvalue"] = _FS.hit_pvalue_get
    if _newclass:pvalue = property(_FS.hit_pvalue_get)
    __swig_getmethods__["evalue"] = _FS.hit_evalue_get
    if _newclass:evalue = property(_FS.hit_evalue_get)
    __swig_getmethods__["zvalue"] = _FS.hit_zvalue_get
    if _newclass:zvalue = property(_FS.hit_zvalue_get)
    __swig_getmethods__["oc_cluster"] = _FS.hit_oc_cluster_get
    if _newclass:oc_cluster = property(_FS.hit_oc_cluster_get)
    __swig_getmethods__["cratio"] = _FS.hit_cratio_get
    if _newclass:cratio = property(_FS.hit_cratio_get)
    __swig_getmethods__["kw_score"] = _FS.hit_kw_score_get
    if _newclass:kw_score = property(_FS.hit_kw_score_get)
    def get_subject(*args): return apply(_FS.hit_get_subject,args)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<C hit instance at %s>" % (self.this,)

class hitPtr(hit):
    def __init__(self,this):
        _swig_setattr(self, hit, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, hit, 'thisown', 0)
        _swig_setattr(self, hit,self.__class__,hit)
_FS.hit_swigregister(hitPtr)

class hit_list(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, hit_list, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, hit_list, name)
    __swig_getmethods__["frag_len"] = _FS.hit_list_frag_len_get
    if _newclass:frag_len = property(_FS.hit_list_frag_len_get)
    __swig_getmethods__["s_db"] = _FS.hit_list_s_db_get
    if _newclass:s_db = property(_FS.hit_list_s_db_get)
    __swig_getmethods__["matrix"] = _FS.hit_list_matrix_get
    if _newclass:matrix = property(_FS.hit_list_matrix_get)
    __swig_getmethods__["range"] = _FS.hit_list_range_get
    if _newclass:range = property(_FS.hit_list_range_get)
    __swig_getmethods__["converted_range"] = _FS.hit_list_converted_range_get
    if _newclass:converted_range = property(_FS.hit_list_converted_range_get)
    __swig_getmethods__["kNN"] = _FS.hit_list_kNN_get
    if _newclass:kNN = property(_FS.hit_list_kNN_get)
    __swig_getmethods__["index_name"] = _FS.hit_list_index_name_get
    if _newclass:index_name = property(_FS.hit_list_index_name_get)
    __swig_getmethods__["alphabet"] = _FS.hit_list_alphabet_get
    if _newclass:alphabet = property(_FS.hit_list_alphabet_get)
    __swig_getmethods__["FS_seqs_total"] = _FS.hit_list_FS_seqs_total_get
    if _newclass:FS_seqs_total = property(_FS.hit_list_FS_seqs_total_get)
    __swig_getmethods__["FS_seqs_visited"] = _FS.hit_list_FS_seqs_visited_get
    if _newclass:FS_seqs_visited = property(_FS.hit_list_FS_seqs_visited_get)
    __swig_getmethods__["FS_seqs_hits"] = _FS.hit_list_FS_seqs_hits_get
    if _newclass:FS_seqs_hits = property(_FS.hit_list_FS_seqs_hits_get)
    __swig_getmethods__["index_seqs_total"] = _FS.hit_list_index_seqs_total_get
    if _newclass:index_seqs_total = property(_FS.hit_list_index_seqs_total_get)
    __swig_getmethods__["seqs_visited"] = _FS.hit_list_seqs_visited_get
    if _newclass:seqs_visited = property(_FS.hit_list_seqs_visited_get)
    __swig_getmethods__["seqs_hits"] = _FS.hit_list_seqs_hits_get
    if _newclass:seqs_hits = property(_FS.hit_list_seqs_hits_get)
    __swig_getmethods__["useqs_visited"] = _FS.hit_list_useqs_visited_get
    if _newclass:useqs_visited = property(_FS.hit_list_useqs_visited_get)
    __swig_getmethods__["useqs_hits"] = _FS.hit_list_useqs_hits_get
    if _newclass:useqs_hits = property(_FS.hit_list_useqs_hits_get)
    __swig_getmethods__["start_time"] = _FS.hit_list_start_time_get
    if _newclass:start_time = property(_FS.hit_list_start_time_get)
    __swig_getmethods__["end_time"] = _FS.hit_list_end_time_get
    if _newclass:end_time = property(_FS.hit_list_end_time_get)
    __swig_getmethods__["search_time"] = _FS.hit_list_search_time_get
    if _newclass:search_time = property(_FS.hit_list_search_time_get)
    __swig_getmethods__["max_hits"] = _FS.hit_list_max_hits_get
    if _newclass:max_hits = property(_FS.hit_list_max_hits_get)
    __swig_getmethods__["actual_seqs_hits"] = _FS.hit_list_actual_seqs_hits_get
    if _newclass:actual_seqs_hits = property(_FS.hit_list_actual_seqs_hits_get)
    __swig_getmethods__["accepted"] = _FS.hit_list_accepted_get
    if _newclass:accepted = property(_FS.hit_list_accepted_get)
    __swig_getmethods__["shape"] = _FS.hit_list_shape_get
    if _newclass:shape = property(_FS.hit_list_shape_get)
    __swig_getmethods__["rate"] = _FS.hit_list_rate_get
    if _newclass:rate = property(_FS.hit_list_rate_get)
    __swig_getmethods__["Zmin"] = _FS.hit_list_Zmin_get
    if _newclass:Zmin = property(_FS.hit_list_Zmin_get)
    def __del__(self, destroy= _FS.delete_hit_list):
        try:
            if self.thisown: destroy(self)
        except: pass
    def print_list(*args): return apply(_FS.hit_list_print_list,args)
    def get_hit(*args): return apply(_FS.hit_list_get_hit,args)
    def get_query(*args): return apply(_FS.hit_list_get_query,args)
    def Z_scores(*args): return apply(_FS.hit_list_Z_scores,args)
    def sort_decr(*args): return apply(_FS.hit_list_sort_decr,args)
    def sort_incr(*args): return apply(_FS.hit_list_sort_incr,args)
    def sort_by_seq(*args): return apply(_FS.hit_list_sort_by_seq,args)
    def sort_by_oc(*args): return apply(_FS.hit_list_sort_by_oc,args)
    def sort_by_evalue(*args): return apply(_FS.hit_list_sort_by_evalue,args)
    def sort_by_cratio(*args): return apply(_FS.hit_list_sort_by_cratio,args)
    def sort_by_kwscore(*args): return apply(_FS.hit_list_sort_by_kwscore,args)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<C hit_list instance at %s>" % (self.this,)

class hit_listPtr(hit_list):
    def __init__(self,this):
        _swig_setattr(self, hit_list, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, hit_list, 'thisown', 0)
        _swig_setattr(self, hit_list,self.__class__,hit_list)
_FS.hit_list_swigregister(hit_listPtr)

class index(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, index, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, index, name)
    __swig_getmethods__["db_name"] = _FS.index_db_name_get
    if _newclass:db_name = property(_FS.index_db_name_get)
    __swig_getmethods__["index_name"] = _FS.index_index_name_get
    if _newclass:index_name = property(_FS.index_index_name_get)
    __swig_getmethods__["s_db"] = _FS.index_s_db_get
    if _newclass:s_db = property(_FS.index_s_db_get)
    __swig_getmethods__["ptable"] = _FS.index_ptable_get
    if _newclass:ptable = property(_FS.index_ptable_get)
    __swig_getmethods__["m"] = _FS.index_m_get
    if _newclass:m = property(_FS.index_m_get)
    __swig_getmethods__["db_no_frags"] = _FS.index_db_no_frags_get
    if _newclass:db_no_frags = property(_FS.index_db_no_frags_get)
    __swig_getmethods__["no_bins"] = _FS.index_no_bins_get
    if _newclass:no_bins = property(_FS.index_no_bins_get)
    __swig_getmethods__["no_seqs"] = _FS.index_no_seqs_get
    if _newclass:no_seqs = property(_FS.index_no_seqs_get)
    __swig_getmethods__["no_useqs"] = _FS.index_no_useqs_get
    if _newclass:no_useqs = property(_FS.index_no_useqs_get)
    def __init__(self,*args):
        _swig_setattr(self, index, 'this', apply(_FS.new_index,args))
        _swig_setattr(self, index, 'thisown', 1)
    def __del__(self, destroy= _FS.delete_index):
        try:
            if self.thisown: destroy(self)
        except: pass
    def save(*args): return apply(_FS.index_save,args)
    def print_stats(*args): return apply(_FS.index_print_stats,args)
    def print_bin(*args): return apply(_FS.index_print_bin,args)
    def get_bin_size(*args): return apply(_FS.index_get_bin_size,args)
    def get_unique_bin_size(*args): return apply(_FS.index_get_unique_bin_size,args)
    def get_seq(*args): return apply(_FS.index_get_seq,args)
    def rng_srch(*args): return apply(_FS.index_rng_srch,args)
    def kNN_srch(*args): return apply(_FS.index_kNN_srch,args)
    def __repr__(self):
        return "<C index instance at %s>" % (self.this,)

class indexPtr(index):
    def __init__(self,this):
        _swig_setattr(self, index, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, index, 'thisown', 0)
        _swig_setattr(self, index,self.__class__,index)
_FS.index_swigregister(indexPtr)

index_load = _FS.index_load

sscan_qd_srch = _FS.sscan_qd_srch

sscan_has_nbr = _FS.sscan_has_nbr


