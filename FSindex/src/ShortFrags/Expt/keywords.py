from Bio.SwissProt import SProt
import sys
import string

class SprotKeywords(dict):
   def  __init__(self):
       dict.__init__(self)
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
