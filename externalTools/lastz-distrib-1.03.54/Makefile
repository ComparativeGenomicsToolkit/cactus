include make-include.mak

default: build_lastz

lastz_32: build_lastz_32

#---------
# builds/installation
#---------

build: build_lastz

build_lastz:
	cd src && ${MAKE} lastz lastz_D

build_lastz_32:
	cd src && ${MAKE} lastz_32

build_test_version:
	cd src && ${MAKE} lastz-test lastz_D-test

clean:
	cd src && ${MAKE} clean

install: install_lastz

install_lastz:
	cd src && ${MAKE} install

install_32:
	cd src && ${MAKE} install_32

install_test_version:
	cd src && ${MAKE} install_test_version

#---------
# testing
#
# Small tests to give some comfort level that the program has built properly,
# or that changes you've made to the source code haven't broken it.  The
# results should be of this form:
#	SUCCESS: ../test_data/xxx and ../test_results/yyy are equivalent
#---------

test:
	cd src && ${MAKE} test

