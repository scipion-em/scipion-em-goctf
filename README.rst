============
goCTF plugin
============

This plugin provides a wrapper for `goCTF <https://www.lsi.umich.edu/science/centers-technologies/cryo-electron-microscopy/research/goctf>`_ program.

.. image:: https://img.shields.io/pypi/v/scipion-em-goctf.svg
        :target: https://pypi.python.org/pypi/scipion-em-goctf
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-goctf.svg
        :target: https://pypi.python.org/pypi/scipion-em-goctf
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-goctf.svg
        :target: https://pypi.python.org/pypi/scipion-em-goctf
        :alt: Supported Python versions

.. image:: https://img.shields.io/sonar/quality_gate/scipion-em_scipion-em-goctf?server=https%3A%2F%2Fsonarcloud.io
        :target: https://sonarcloud.io/dashboard?id=scipion-em_scipion-em-goctf
        :alt: SonarCloud quality gate

.. image:: https://img.shields.io/pypi/dm/scipion-em-goctf
        :target: https://pypi.python.org/pypi/scipion-em-goctf
        :alt: Downloads

Installation
------------

You will need to use 3.0+ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

    scipion installp -p scipion-em-goctf

b) Developer's version

    * download repository

    .. code-block::

        git clone https://github.com/scipion-em/scipion-em-goctf.git

    * install

    .. code-block::

        scipion installp -p /path/to/scipion-em-goctf --devel

goCTF binary will be installed automatically with the plugin, but you can also link an existing installation.
Default installation path assumed is ``software/em/goctf-1.2.0``, if you want to change it, set *GOCTF_HOME* in ``scipion.conf`` file to
the folder where the goCTF is installed.

To check the installation, simply run the following Scipion test:

``scipion tests goctf.tests.test_protocols_goctf.TestGoCTF``


Supported versions
------------------

1.2.0

Protocols
---------

    * ctf refinement


References
----------

    * Min Su (2019). goCTF: Geometrically optimized CTF determination for single-particle cryo-EM. JSB 205(1), 22-29. https://doi.org/10.1016/j.jsb.2018.11.012
