.. the underscore + name + ':' creates a flag for referring to this RST file in other RST files/this file

.. _VHVLPairing:

.. here is a title
 
VH-VL Pairing
=============

The following describes the current methods we use for analyzing our VH-VL Paired sequence data.

.. formatting headers -> equal signs refers to HTML <H1> tag whereas apostrophes ''' refers to <h3> tag
.. BULLET POINTS -> USE ASTERICES
.. using references -> :ref: label for reference <its actual link to a flag> 
.. using automodules => members gives all document function. So function with """ and """ (docstrings);
		-> undoc-members will just list undocumentated functions 

Contents
''''''''
* :ref:`Pairing functions <pairingfunc>`
* :ref:`Single VH chain clustering <singlechain>`

* Examples for usage


.. _pairingfunc:

.. automodule:: immunogrep_gglab_pairing
	:members:
	:undoc-members: populate_clusters, gene_hist, ProcessGene,WriteSummaryFile

	.. _singlechain:
	
	.. note:: This analysis is a supplement to VH-VL pairing. It performs the same analysis as VH-VL clustering excepts omits light chain data.

.. automodule:: immunogrep_singlechain_cdr3_cluster
	:members:
