import FS

class db(object):
    def __init__(self, db_name, new=True, own=True):
        if new == True:
            self.this = FS.new_db(db_name)
            self.thisown = 1
        elif own == True:
            # SWIG string - thisown = 1
            self.this = db_name
            self.thisown = 1
        else:
            # SWIG string - thisown = 0
            self.this = db_name
            self.thisown = 0
            
    def __del__(self):
        if self.thisown:
            FS.delete_db(self) 

    def __str__(self):
      return 'FASTA SEQUENCE DATABASE\nFile: %s\nLength: %d\nSequences: %d\n' % (self.db_name, self.length , self.no_seq)    

    def get_seq(self, i):
        return FS.db_get_seq(self, i)

    def get_def(self, i):
        return FS.db_get_def(self, i)

    def get_frag(self, i, a, b):
        return FS.db_get_frag(self, i, a, b)

    db_name = property(FS.db_db_name_get)
    length = property(FS.db_length_get)
    no_seq = property(FS.db_no_seq_get)
