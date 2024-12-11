Aircraft Geometry Generation
============================

Overview
--------

AGG is a library and associated tools for generating aircraft
geometries, using the CST methods of Brenda Kulfan. The intention is
to allow users to supply geometries in a parametric form which allows
automatic generation of surface geometries which can be used in
computations.

After installing prerequisites, installation follows the standard
procedure:

./configure && make && make install

Prerequisites
-------------

You will need to install the following GNU tools:

    autoconf
    automake

AGG now requires the libmatheval library for evaluation and
differentiation of symbolic expressions. This is available from

https://www.gnu.org/software/libmatheval/

but the test suite for this version requires (an old version of) Guile
which may not be available on some systems. A fork of the project
which does not require Guile is available from

https://github.com/mjcarley/libmatheval-no-guile

AGG uses Jonathan Shewchuk's Triangle code for mesh generation,
through Christian Woltering's API, available from:

https://github.com/wo80/Triangle/

You should install this before installing AGG.

Step-by-Step
------------

Clone this repository:

    git clone https://github.com/mjcarley/agg.git

Generate the configure script:

    autoreconf -ivf

Configure and build the project:

    ./configure
    make
