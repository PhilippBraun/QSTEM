# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/cmake-gui

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/philipp/QSTEM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/philipp/QSTEM

# Include any dependencies generated for this target.
include tests/libs/CMakeFiles/test_detectors.dir/depend.make

# Include the progress variables for this target.
include tests/libs/CMakeFiles/test_detectors.dir/progress.make

# Include the compile flags for this target's objects.
include tests/libs/CMakeFiles/test_detectors.dir/flags.make

tests/libs/CMakeFiles/test_detectors.dir/test_detectors.cpp.o: tests/libs/CMakeFiles/test_detectors.dir/flags.make
tests/libs/CMakeFiles/test_detectors.dir/test_detectors.cpp.o: tests/libs/test_detectors.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/philipp/QSTEM/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object tests/libs/CMakeFiles/test_detectors.dir/test_detectors.cpp.o"
	cd /home/philipp/QSTEM/tests/libs && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test_detectors.dir/test_detectors.cpp.o -c /home/philipp/QSTEM/tests/libs/test_detectors.cpp

tests/libs/CMakeFiles/test_detectors.dir/test_detectors.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_detectors.dir/test_detectors.cpp.i"
	cd /home/philipp/QSTEM/tests/libs && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/philipp/QSTEM/tests/libs/test_detectors.cpp > CMakeFiles/test_detectors.dir/test_detectors.cpp.i

tests/libs/CMakeFiles/test_detectors.dir/test_detectors.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_detectors.dir/test_detectors.cpp.s"
	cd /home/philipp/QSTEM/tests/libs && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/philipp/QSTEM/tests/libs/test_detectors.cpp -o CMakeFiles/test_detectors.dir/test_detectors.cpp.s

tests/libs/CMakeFiles/test_detectors.dir/test_detectors.cpp.o.requires:
.PHONY : tests/libs/CMakeFiles/test_detectors.dir/test_detectors.cpp.o.requires

tests/libs/CMakeFiles/test_detectors.dir/test_detectors.cpp.o.provides: tests/libs/CMakeFiles/test_detectors.dir/test_detectors.cpp.o.requires
	$(MAKE) -f tests/libs/CMakeFiles/test_detectors.dir/build.make tests/libs/CMakeFiles/test_detectors.dir/test_detectors.cpp.o.provides.build
.PHONY : tests/libs/CMakeFiles/test_detectors.dir/test_detectors.cpp.o.provides

tests/libs/CMakeFiles/test_detectors.dir/test_detectors.cpp.o.provides.build: tests/libs/CMakeFiles/test_detectors.dir/test_detectors.cpp.o

# Object files for target test_detectors
test_detectors_OBJECTS = \
"CMakeFiles/test_detectors.dir/test_detectors.cpp.o"

# External object files for target test_detectors
test_detectors_EXTERNAL_OBJECTS =

bin/test_detectors: tests/libs/CMakeFiles/test_detectors.dir/test_detectors.cpp.o
bin/test_detectors: tests/libs/CMakeFiles/test_detectors.dir/build.make
bin/test_detectors: bin/libqstem.so
bin/test_detectors: /usr/local/lib/libfftw3.so
bin/test_detectors: /usr/lib/x86_64-linux-gnu/libfftw3f.so
bin/test_detectors: /usr/local/lib/libboost_unit_test_framework.so
bin/test_detectors: /usr/local/lib/libboost_system.so
bin/test_detectors: /usr/local/lib/libboost_filesystem.so
bin/test_detectors: /usr/local/lib/libhdf5.so
bin/test_detectors: /usr/lib/x86_64-linux-gnu/libpthread.so
bin/test_detectors: /usr/lib/x86_64-linux-gnu/libdl.so
bin/test_detectors: /usr/lib/x86_64-linux-gnu/libm.so
bin/test_detectors: /usr/local/lib/libfftw3.so
bin/test_detectors: /usr/lib/x86_64-linux-gnu/libfftw3f.so
bin/test_detectors: /usr/local/lib/libfftw3_omp.a
bin/test_detectors: /usr/local/lib/libcvmlcpp.so
bin/test_detectors: /usr/local/lib/libblitz.a
bin/test_detectors: tests/libs/CMakeFiles/test_detectors.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/test_detectors"
	cd /home/philipp/QSTEM/tests/libs && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_detectors.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/libs/CMakeFiles/test_detectors.dir/build: bin/test_detectors
.PHONY : tests/libs/CMakeFiles/test_detectors.dir/build

tests/libs/CMakeFiles/test_detectors.dir/requires: tests/libs/CMakeFiles/test_detectors.dir/test_detectors.cpp.o.requires
.PHONY : tests/libs/CMakeFiles/test_detectors.dir/requires

tests/libs/CMakeFiles/test_detectors.dir/clean:
	cd /home/philipp/QSTEM/tests/libs && $(CMAKE_COMMAND) -P CMakeFiles/test_detectors.dir/cmake_clean.cmake
.PHONY : tests/libs/CMakeFiles/test_detectors.dir/clean

tests/libs/CMakeFiles/test_detectors.dir/depend:
	cd /home/philipp/QSTEM && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/philipp/QSTEM /home/philipp/QSTEM/tests/libs /home/philipp/QSTEM /home/philipp/QSTEM/tests/libs /home/philipp/QSTEM/tests/libs/CMakeFiles/test_detectors.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/libs/CMakeFiles/test_detectors.dir/depend

