Running FSIndex search servers
==============================

There are two types of an FSIndex server: a slave, which loads an index and conducts searches and a master, which distributes and collects searches. They communicate between themselves and with clients using TCP/IP sockets.

The master/slave configuration can be used in the case where indexes are located on different machines, even running different operating systems. If there is only one server, it runs under the slave configuration. This means that each slave can be queried individually and also that a more complicated tree-like structure, where inner nodes are master servers and leaves are slaves, can be constructed if needed.

Two scripts can be used to run an FSIndex server: ``FSsearchd.py``, and ``FSsearchc.py``. ``FSsearchd.py`` runs as a UNIX daemon, while ``FSsearchc.py`` is a regular script for all platforms. Their arguments and options are similar and are described in detail only for ``FSsearchd.py``.

.. _subsec-FSsearchd.py:

FSsearchd.py daemon
-------------------

FSsearchd.py script starts an FSIndex server as a UNIX daemon. Unlike the other script, the master slave can be passed an option to control its slaves, in which case it starts the slaves using ``ssh`` on startup and shut them down when it terminates.. The current version assumes that the ``ssh`` login is using public key without a password.

Starting
^^^^^^^^

``FSsearchd.py`` daemon is started from the command line::

    $ FSsearchd.py [-c] serverid port workpath indexfile [pythonpath]

where

* ``serverid`` is a string identifier of a particular index (or set of indexes) to be loaded;

* ``port`` is the port the daemon should be listening to;

* ``workpath`` is the path to the directory where log files are written;

* ``indexfile`` is the path to the index file. If the file ends with ``.cfg``, it is assumed to be a configuration file for the slaves;

* ``pythonpath`` (optional) is appended to the system's python path. This may help in the case that the path to Python modules required by ``FSsearchd.py`` is not set otherwise.

The option ``--control-slaves`` or ``-c`` can be used to instruct a master server to attempt to start and terminate its slaves using ssh. If this option is omitted, the master will only attempt to contact each slave and, if any one is unavailable, shut itself down. In this case the slaves must be started separately.

A master's configuration file is a text file where each line gives parameters for a slave server to be started. The format of each line is::

    host port workpath indexfile [pythonpath [binpath]]

where the fields are separated by spaces. The fields ``port``, ``workpath``, ``indexfile`` and ``pythonpath`` are used directly to start the slave server (see above for full description). The field ``host`` specifies the address of the machine the slave is running on. The optional parameter ``binpath`` provides the path to the ``FSsearchd.py`` executable on the host the slave should run on. Only ``host`` and ``port`` fields are used if ``-c`` option is not set.

Each server produces two files: a log file and a pid file. The log file receives detailed messages about the running of the server while the pid file contains the UNIX process id of the daemon, the name of the host it is running on, as well as the command line used to start it. The name of the log file is generated from the ``serverid`` passed as a command line argument. If a master is starting its slaves, it passes its own ``serverid`` to them, concatenated with the number of the slave. For example, a master server passed ``TEST2`` as ``serverid`` will produce the logfile named ``FSsearchd_TEST2.log`` and its second slave will produce ``FSsearchd\_TEST2\_s01.log``.

Terminating
^^^^^^^^^^^

Since ``FSsearchd.py`` is a daemon and hence not connected to any terminal, the best way to terminate it is to send it a SIGTERM signal. To do so, find out (from a pid file or the output of \textit{ps}) its process id (``pid`` and type::

    $ kill pid

This is a clean way to shutdown an FSIndex server since the logs are written and any pending requests are handled before shutdown.


FSsearchc.py script
-------------------

``FSsearchc.py`` script runs an FSIndex server from the command line, printing log to the standard output. It should run on all platforms and is started by::

    $ FSsearchc.py port indexfile

where the parameters are as described in :ref:`subsec-FSsearchd.py` above. It can be stopped by killing the server process.
