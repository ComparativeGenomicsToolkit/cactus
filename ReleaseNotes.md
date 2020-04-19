
# Release 1.0.0   2020-04-19

This is the first official release of Cactus.  The goal is to provide a
stable, track-able version of Cactus to users.  The releases is provided in
source, static-compiled binaries, and Docker formats.

Notable Changes:
 - Kyoto Cabinet and Typhoon are now included as a sub-module.  This is due to
   the lack of consistent, stable releases and the difficulty in compiling it.
 - Cactus is now available as a tar file of static-linked binary executables,
   along with a wheel of the Cactus Python packages.  This avoids compilation and dependency problems. These should work on most Linux platforms when Docker is not used.
 - Added support to run Cactus in individual steps. This works around problems with using Toil in some distributed environments by dividing up alignment into tasks that can be run manually on separate machines.
   See *Running step by step (experimental)* in `README.md`.
 - Conversion to Python 3, allowing Toil to drop Python 2 support.


