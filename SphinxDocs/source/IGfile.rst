
.. _immunogrepfile:

Immunogrepfile
==============

Contents:

* :ref:`Variables <filevariables>`
* :ref:`Header Extraction Functions <fileheaders>`
    * :py:func:`.GetAdditionalInfo`
    * :py:func:`.find_imgt_file_type_index`
    * :py:func:`.GroupIMGTFiles`
    * :py:func:`.ReadIgBlastQueryBlock`
* :ref:`Classes <fileclasses>`
    * :py:class:`.immunogrepFile`
    * :py:class:`.immunogrepFASTA`
    * :py:class:`.immunogrepFASTQ`
    * :py:class:`.immunogrepDELIM`
    * :py:class:`.immunogrepIMGT`
    * :py:class:`.immunogrepJSON`

.. _filevariables:

Variables
'''''''''
:current_filetypes:

    =========   =================
    Filetype
    ---------   -----------------
    FASTA       fna, fasta
    TAB         txt, tab
    CSV         csv
    FASTQ       fastq
    JSON        json, annotation, analysis, query, iffile, IFFile
    DELIM       delim
    IMGT        imgt
    =========   =================

:ref:`^To Top <immunogrepfile>`

.. _fileheaders:

Header Extraction Functions
'''''''''''''''''''''''''''

.. automodule:: immunogrep_read_file
    :members: GetAdditionalInfo, find_imgt_file_type_index, GroupIMGTFiles, ReadIgBlastQueryBlock
    :undoc-members:

:ref:`^To Top <immunogrepfile>`

.. _fileclasses:

Immunogrep Classes
''''''''''''''''''
    Consists of several classes to import various file types. Many of the classes have similar functions.

    .. note::
        :py:class:`.immunogrepFile` is the primary class that is run.

.. autoclass:: immunogrepFile
    :members:
    :undoc-members:

.. autoclass:: immunogrepFASTA
    :members:
    :undoc-members:

.. autoclass:: immunogrepFASTQ
    :members:
    :undoc-members:

.. autoclass:: immunogrepDELIM
    :members:
    :undoc-members:

.. autoclass:: immunogrepJSON
    :members:
    :undoc-members:

.. autoclass:: immunogrepIMGT
    :members:
    :undoc-members:

:ref:`^To Top <immunogrepfile>`



