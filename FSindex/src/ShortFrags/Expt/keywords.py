from Bio.SwissProt import SProt
import sys, string, shelve
from cStringIO import StringIO

class CachedShelve(shelve.DbfilenameShelf):
   def __init__(self, *args, **kwargs):
      shelve.DbfilenameShelf(self, *args, **kwargs)
      self._cache = {}
    
   def __getitem__(self, key):
      return self._cache.setdefault(key, self[key])

   def __delitem__(self, key):
      if key in self._cache:
         del(self._cache[key])
      del(self[key])

   def __contains__(self, key):
      if key in self._cache or key in self:
         return True
      else:
         return False

   def has_value(self, value):
      if value in self._cache.values():
         pass

   def cache_len(self):
      return len(self._cache)

   def clear_cache(self):
      self._cache = {}


## class SequenceKeywords(object):
##    def __init__(self, *args, **kwargs):
##       self._map = CachedShelve()
##       self._values = CachedShelve()
##       if flag != 'r':
##          self._no_insert = False
##          val_ids = self._values.keys()
##          self._counts =  
##          pass
##       else:
##          self._no_insert = True
            
##       self._tot_vals = len(self._values)

         
##    def __getitem__(self, key):
##       vals0 = self._map[key]
##       value = []
##       for val_id in vals0:
##          value.append[self._values[val_id][0]]
##       return value

##    def __setitem__(self, key, value):
##       if self._no_insert:
##          throw RuntimeError: "No insertion allowed."
##       for val in value:
##          if val in self._val_dict:
##             val_id = self._val_dict[val]
##             (v,c) = self._values[val_id]
##             c += 1
##             self._values[val_id] = (v,c)
##          else:
##             val_id = self._tot_vals + 1
##             self._val_dict[val] = val_id
##             self._values[val_id] = (value, 1)
##             self._tot_vals += 1


       
##          def __delitem__(self, key):
##         if key in self._cache:
##             del(self._cache[key])
##         del(self[key])

##     def cache_len(self):
##         return len(self._cache)

##     def clear_cache(self):
##         self._cache = {}



class SprotKeywords(dict):
   def  __init__(self):
      self.kw_dict = {}
      self.kw = list()
      self.freq = list() 

   def parse_file(self, filename):
      fp = file(filename, 'r')
      s_parser = SProt.RecordParser()
      s_iterator = SProt.Iterator(fp, s_parser)

      i = 0
      kwid = 0
      for rec in s_iterator:
         kwl = list()
         for k in rec.keywords:
            if k in self.kw_dict:
               self.kw_dict[k][1] += 1
               kwl.append(self.kw_dict[k][0])
            else:
               self.kw_dict[k] = [kwid, 1]
               kwl.append(kwid)
               kwid += 1
         if len(kwl) > 0:
            for ac in rec.accessions:
               self[ac] = kwl
         i += 1
         if i % 1000 == 0:
            print i
      fp.close()      
      self.kw = list(kwid*[0])
      self.freq = list(kwid*[0])
      for k in self.kw_dict:
         kwid = self.kw_dict[k][0]
         self.kw[kwid] = k;
         self.freq[kwid] = self.kw_dict[k][1]
      self.kw_dict = {}
       
   def save(self, filename):
      fp = file(filename, 'w')
      fp.write('# Numbers: keywords, accessions\n')
      fp.write('%d %d\n' % (len(self.kw), len(self)))
      fp.write('# Keywords: id, frequency, keyword\n')
      for (i, k) in enumerate(self.kw):
         fp.write('%d %d %s\n' % (i, self.freq[i], k))
      fp.write('# Accessions: accession, no_keywords, list(keyword_id)\n')

      ackeys = self.keys()
      ackeys.sort()
      for ac in ackeys:
         fp.write('%s %d' % (ac, len(self[ac])))
         for j in self[ac]:
            fp.write(' %d' % j)
         fp.write('\n')
      fp.close()

   def load(self, filename):
      fp = file(filename, 'r')

      fp.next() # comment
      fp.next() # Numbers: not needed
      fp.next() # comment

      while 1: # keywords
         line = fp.next()
         if line[0] == '#':
            break # comment for acc section
         s = string.split(line, None, 2)
         self.freq.append(int(s[1]))
         self.kw.append(s[2][:-1])
           
      for line in fp:
         s = string.split(line)
         self[s[0]] = [int(kwids) for kwids in s[2:]]
      fp.close()

   def get_keywords(self, acc):
      if acc in self:
         return [self.kw[j] for j in self[acc]]
      else:
         return []

class UnirefClusters(SprotKeywords):
   def __init__(self):
      SprotKeywords.__init__(self)
      self._name_flag = False
      self._name = ''
      self._kwid = 0
      
   def _start_element(self, name, attrs):
      if name == 'name':
         self._name_flag = True
         self._fs = StringIO()
         self._kwl = []
      elif name == 'property' and attrs['type'] == "UniProt accession":
         ac = str(attrs['value'])
         self[ac] = self._kwl

   def _end_element(self, name):
      if name == 'name':
         self._name_flag = False
         k = self._fs.getvalue()
         if k in self.kw_dict:
            self.kw_dict[k][1] += 1
            self._kwl.append(self.kw_dict[k][0])
         else:
            self.kw_dict[k] = [self._kwid, 1]
            self._kwl.append(self._kwid)
            self._kwid += 1
##          if self._kwid < 50:
##             print k

   def _char_data(self, data):
      if self._name_flag:
         self._fs.write(str(data))
         
   def parse_file(self, filename):
      import xml.parsers.expat
      p = xml.parsers.expat.ParserCreate()

      p.StartElementHandler = self._start_element
      p.EndElementHandler = self._end_element
      p.CharacterDataHandler = self._char_data


      fp = file(filename, 'r')
      p.ParseFile(fp) 
      fp.close()

      self.kw = list(self._kwid*[0])
      self.freq = list(self._kwid*[0])
      for k in self.kw_dict:
         kwid = self.kw_dict[k][0]
         self.kw[kwid] = k;
         self.freq[kwid] = self.kw_dict[k][1]
      self.kw_dict = {}
