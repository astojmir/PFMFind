import os.path, xml.parsers.expat
from cStringIO import StringIO

from BioSQL import BioSeqDatabase

from ShortFrags.Setup.UniprotLoader import load_Uniprot
from ShortFrags.Setup.UnirefLoader import UnirefParser
from ShortFrags.Setup.UnirefLoader import UnirefLoader
from ShortFrags.Setup.InterProLoader import Protein2InterProParser
from ShortFrags.Setup.InterProLoader import InterProLoader

class PFMF_DatabaseLoader(object):
    def __init__(self):
        self._setup_flag = False
        self._sflag = False
        
        self.dbargs = {}
        self.schema_name = None
        self.schema_create = 0
        self.initial_sql = None
        self.final_sql = None

        self.download_taxonomy = 0
        self.copy_taxonomy = 0
        self.perl_driver = None
        self.taxon_dir = None
        self.script_dir = './'

        self.uniprot = []
        self.uniref = []
        self.interpro = []

    def _start_element(self, name, attrs):

        if name == 'PFMF_db_setup':
            self._setup_flag = True 

        if not self._setup_flag: return
        
        if name == 'Database':
            self.dbargs = attrs
        elif name == 'Uniprot_file' or name == 'Uniref_file' \
                 or name == 'InterPro_file':
            self._namespace = attrs['namespace']
            self._sflag = True
            self._fs = StringIO()
            self._sql_start = attrs.get('sql_start', None)
            self._sql_end = attrs.get('sql_end', None)
        elif name == 'Sql_dir' or name == 'Taxon_dir' \
                 or name == 'Script_dir':
            self._sflag = True
            self._fs = StringIO()
        elif name == 'Sql_scripts':
            self.initial_sql = attrs.get('sql_start', None)
            self.final_sql = attrs.get('sql_end', None)
        elif name == 'Schema':
            self.schema_name = attrs.get('name', None)
            self.schema_create = int(attrs.get('create', '0'))
        elif name == 'Taxonomy':
            self.download_taxonomy = int(attrs.get('download', '0'))
            self.perl_driver = attrs.get('perl_driver', None)
            self.copy_taxonomy = int(attrs.get('copy', '0'))
            
    def _end_element(self, name):

        if not self._setup_flag: return

        if name == 'PFMF_db_setup':
            self._setup_flag = False 
        elif name == 'Uniprot_file':
            self._sflag = False
            self.uniprot.append((self._fs.getvalue(), self._namespace,
                                 self._sql_start, self._sql_end)) 
        elif name == 'Uniref_file':
            self._sflag = False
            self.uniref.append((self._fs.getvalue(), self._namespace, 
                                self._sql_start, self._sql_end))
        elif name == 'InterPro_file':
            self._sflag = False
            self.interpro.append((self._fs.getvalue(),
                                  self._namespace,
                                  self._sql_start, self._sql_end))
        elif name == 'Sql_dir':
            self._sflag = False
            self.sql_dir = self._fs.getvalue()
        elif name == 'Taxon_dir':
            self._sflag = False
            self.taxon_dir = self._fs.getvalue()
        elif name == 'Script_dir':
            self._sflag = False
            self.script_dir = self._fs.getvalue()

    def _char_data(self, data):

        if self._setup_flag and self._sflag:
            self._fs.write(data)
         
    def parse_config(self, fp):
        """
        Parses the XML configuration file.
        """

        p = xml.parsers.expat.ParserCreate()
        p.returns_unicode = 0
        p.StartElementHandler = self._start_element
        p.EndElementHandler = self._end_element
        p.CharacterDataHandler = self._char_data
        
        p.ParseFile(fp) 

    def load_database(self):
        """
        Loads a relational database according to previously
        parsed configuration.
        """

        # Open database
        print "Opening database."
        self._init_database()

        # Setting Schema
        if self.schema_name:
            print "Setting up the schema %s." % self.schema_name 
            self._set_schema()

        # Initial SQL
        if self.initial_sql:
            print "Running initial SQL script %s." % self.initial_sql
            self._execute_sql(self.initial_sql)
            
        # Taxonomy
        if self.taxon_dir:
            print "Loading taxonomy."
            if self.copy_taxonomy:
                self._copy_taxonomy()
            else:
                self._load_taxonomy()
            
        # Uniprot
        for filename, namespace, sql0, sql1 in self.uniprot:
            print "Processing Uniprot file %s"\
                  % os.path.split(filename)[1]
            if sql0:
                print "    Running initial SQL script %s." % sql0
                self._execute_sql(sql0)
            print "    Loading .dat file ..."
            n = self._load_uniprot(filename, namespace)
            print "    Loaded %d records." % n
            if sql1:
                print "    Running final SQL script %s." % sql1
                self._execute_sql(sql1)

        # Uniref
        for filename, namespace, sql0, sql1 in self.uniref:
            print "Processing Uniref file %s"\
                  % os.path.split(filename)[1]
            if sql0:
                print "    Running initial SQL script %s." % sql0
                self._execute_sql(sql0)
            print "    Loading XML file ..."
            stats = self._load_uniref(filename, namespace)
            print "    Parsed %d, %d. Loaded %d, %d." % stats
            if sql1:
                print "    Running final SQL script %s." % sql1
                self._execute_sql(sql1)

        # InterPro
        for filename, namespace, sql0, sql1 in self.interpro:
            print "Processing InterPro file %s"\
                  % os.path.split(filename)[1]
            if sql0:
                print "    Running initial SQL script %s." % sql0
                self._execute_sql(sql0)
            print "    Loading protein2interpro file ..."
            stats = self._load_interpro(filename, namespace)
            print "    Lines %s, Domains %s." % stats
            if sql1:
                print "    Running final SQL script %s." % sql1
                self._execute_sql(sql1)

        # Final SQL
        if self.final_sql:
            print "Running final SQL script %s." % self.final_sql
            self._execute_sql(self.final_sql)

        # Grant reading priviledges on schema
        if self.schema_name:
            print "Granting reading priviledges on schema."
            self._make_schema_readable()

        print "Closing database."
        self.server.adaptor.close()
        
    def _execute_sql(self, sql_file):
        cur = self.server.adaptor.cursor
        fn1 = os.path.join(self.sql_dir, sql_file)
        fp = file(fn1, 'r')
        s = fp.read()
        fp.close()
        cur.execute(s)
        self.server.adaptor.commit()

    def _init_database(self):
        """
        Opens a relational database.
        """
        
        self.server = BioSeqDatabase.open_database(**self.dbargs)


    def _set_schema(self):
        """
        Sets current schema (optionally creates) - PostgreSQL only.
        """

        cur = self.server.adaptor.cursor
        if self.schema_create:
            cur.execute("CREATE SCHEMA %s" % self.schema_name)

        cur.execute("SET search_path TO %s" % self.schema_name)
        self.server.adaptor.commit()

    def _make_schema_readable(self):
        """
        Makes schema readable for all users - PostgreSQL only.
        """
        
        cur = self.server.adaptor.cursor
        cur.execute('GRANT USAGE ON SCHEMA %s TO PUBLIC'\
                    % (self.schema_name,)) 
        sql = "SELECT table_name FROM information_schema.tables" \
              " WHERE table_schema=%s"  
        cur.execute(sql, (self.schema_name,))

        cur1 = self.server.adaptor.conn.cursor()
        while 1:
            tn = cur.fetchone()
            if tn == None: break
            cur1.execute('GRANT SELECT ON TABLE %s TO PUBLIC' % tn[0])
        cur1.close()

        self.server.adaptor.commit()

    def _copy_taxonomy(self):
        """
        Loads taxonomy tables from tab-separeted files (PostgreSQL only).
        """
        cur = self.server.adaptor.cursor
        taxon_tbl = os.path.join(self.taxon_dir, 'taxon.tbl')
        taxon_name_tbl = os.path.join(self.taxon_dir, 'taxon_name.tbl')
        
        cur.execute("COPY taxon FROM %s", (taxon_tbl,))
        cur.execute("COPY taxon_name FROM %s", (taxon_name_tbl,))
        self.server.adaptor.commit()

    def _load_taxonomy(self):
        """
        Loads taxonomy tables using the load_ncbi_taxonomy.pl script.
        """

        # Set search path to schema - PostgreSQL only
        if self.schema_name != None:
            cur = self.server.adaptor.cursor
            cur.execute('SHOW search_path')
            old_db_path = cur.fetchone()[0]
            cur.execute('SELECT SESSION_USER')
            user = cur.fetchone()[0]

            cur.execute('ALTER USER %s SET search_path TO %s' \
                        % (user, self.schema_name))
            
            self.server.adaptor.commit()

        script = os.path.join(self.script_dir, 'load_ncbi_taxonomy.pl')
        args = [script,
                '--driver=%s' % self.perl_driver,
                '--dbname=%s' % self.dbargs['db'],
                '--verbose=2']

        if 'host' in self.dbargs:
            args.append('--host=%s' % self.dbargs['host'])
        if 'port' in self.dbargs:
            args.append('--port=%s' % self.dbargs['port'])
        if 'user' in self.dbargs:
            args.append('--dbuser=%s' % self.dbargs['user'])
        if 'password' in self.dbargs:
            args.append('--dbpass=%s' % self.dbargs['password'])
        if self.taxon_dir:
            args.append('--directory=%s' % self.taxon_dir)
        if self.download_taxonomy:
            args.append('--download')
    
        os.spawnv(os.P_WAIT, script, args)
    
        # Reset the search path
        if self.schema_name != None:
            cur.execute('ALTER USER %s SET search_path TO %s',
                        (user, old_db_path))

            self.server.adaptor.commit()


    def _get_namespace(self, namespace):
        if namespace in self.server.keys():
            return self.server[namespace]
        else:
            return self.server.new_database(namespace)

    def _load_uniprot(self, filename, namespace):
        """
        Loads a Uniprot .dat file.
        """

        nsdb = self._get_namespace(namespace)

        n = load_Uniprot(nsdb.adaptor, nsdb.dbid, filename)
        self.server.adaptor.commit()
        return n

    def _load_uniref(self, filename, namespace):
        """
        Loads Uniref clusters from XML file.
        """

        nsdb = self._get_namespace(namespace)
        uniref_name = os.path.split(filename)[1]

        UL = UnirefLoader(nsdb.adaptor, nsdb.dbid, uniref_name)
        UP = UnirefParser(UL.load_cluster)
        fp = file(filename)
        UP.parse(fp)
        fp.close()
        self.server.adaptor.commit()
        return (UP.total_clusters, UP.total_accessions,
                UL.total_clusters, UL.total_accessions)

    def _load_interpro(self, filename, namespace):
        """
        Loads InterPro domains from protein2interpro.dat file.
        """

        nsdb = self._get_namespace(namespace)

        IL = InterProLoader(nsdb.adaptor, nsdb.dbid)
        IP = Protein2InterProParser(IL.load_domain)
        fp = file(filename)
        IP.parse(fp)
        fp.close()
        self.server.adaptor.commit()
        return (IP.total_lines, IL.total_domains)
